"""Passive, globally weighted monthly RRTMGP cloud-layer summaries."""
mutable struct CloudProfileCapture
    start::DateTime
    stop::DateTime
    weights::Vector{Float64}
    sigma::Vector{Float64}
    dsigma::Vector{Float64}
    scheme::String
    fraction_sum::Vector{Float64}
    liquid_sum::Vector{Float64}
    ice_sum::Vector{Float64}
    top_sum::Vector{Float64}
    samples::Int
end

function CloudProfileCapture(start, stop, weights, sigma, dsigma, scheme)
    nlayer = length(sigma)
    nlayer > 0 && length(dsigma) == nlayer || error("invalid cloud layer coordinates")
    all(isfinite, weights) && all(>=(0), weights) &&
        isapprox(sum(weights), 1; atol=1e-5, rtol=0) ||
        error("cloud profile point weights do not close")
    all(isfinite, sigma) && all(isfinite, dsigma) && all(>(0), dsigma) ||
        error("invalid cloud sigma coordinates")
    normalized_weights = Float64.(weights) ./ sum(weights)
    return CloudProfileCapture(start, stop, normalized_weights, Float64.(sigma), Float64.(dsigma),
        scheme, zeros(nlayer), zeros(nlayer), zeros(nlayer), zeros(nlayer + 1), 0)
end

function _cloud_input_hashes(vars, model)
    hashes = R._scored_input_hashes(vars)
    radiation = model.longwave_radiation
    if radiation isa R.RRTMGPRadiation && radiation.cloud_scheme in R._ALL_SKY_CLOUD_SCHEMES
        solver = radiation.solver
        for values in (R.RRTMGP.cloud_fraction(solver),
                R.RRTMGP.cloud_liquid_water_path(solver),
                R.RRTMGP.cloud_ice_water_path(solver))
            push!(hashes, bytes2hex(sha256(reinterpret(UInt8, vec(Array(values))))))
        end
    end
    return hashes
end

function _cloud_sample!(capture::CloudProfileCapture, vars, model, time::DateTime)
    capture.start <= time < capture.stop || error("cloud state sample outside window")
    time == capture.start + Hour(capture.samples) ||
        error("cloud state sample missing, duplicated or not hourly")
    before = _cloud_input_hashes(vars, model)
    fields = R._atmosphere_cloud_layer_diagnostics(model)
    fraction, liquid, ice = fields.fraction, fields.liquid_water_path,
        fields.ice_water_path
    npoint, nlayer = size(fraction)
    npoint == length(capture.weights) && nlayer == length(capture.sigma) &&
        size(liquid) == size(ice) == size(fraction) &&
        length(fields.column_fraction) == npoint || error("cloud profile shape differs")
    all(isfinite, fraction) && all(x -> 0 <= x <= 1, fraction) &&
        all(isfinite, liquid) && all(>=(0), liquid) &&
        all(isfinite, ice) && all(>=(0), ice) &&
        all(isfinite, fields.column_fraction) &&
        all(x -> 0 <= x <= 1, fields.column_fraction) ||
        error("cloud profile contains invalid values")
    sample_bins = zeros(nlayer + 1)
    for point in 1:npoint
        weight = capture.weights[point]
        cover = fields.column_fraction[point]
        top = findfirst(>(0), view(fraction, point, :))
        isnothing(top) ? (cover == 0 || error("cloud cover without a cloud layer")) :
            (sample_bins[top] += weight * cover)
        sample_bins[end] += weight * (1 - cover)
        for layer in 1:nlayer
            cloudy = fraction[point, layer]
            capture.fraction_sum[layer] += weight * cloudy
            capture.liquid_sum[layer] += weight * cloudy * liquid[point, layer]
            capture.ice_sum[layer] += weight * cloudy * ice[point, layer]
        end
    end
    all(isfinite, sample_bins) && all(>=(0), sample_bins) &&
        isapprox(sum(sample_bins), 1; atol=1e-12, rtol=0) ||
        error("cloud-top histogram does not close")
    capture.top_sum .+= sample_bins
    capture.samples += 1
    before == _cloud_input_hashes(vars, model) ||
        error("cloud profile sampling changed model inputs")
    return nothing
end

function _cloud_finish(capture::CloudProfileCapture)
    expected = Dates.value(capture.stop - capture.start) ÷ 3_600_000
    capture.samples == expected || error("cloud profile window has incomplete hourly coverage")
    fraction = capture.fraction_sum ./ expected
    liquid = capture.liquid_sum ./ expected
    ice = capture.ice_sum ./ expected
    top = capture.top_sum ./ expected
    all(isfinite, fraction) && all(x -> 0 <= x <= 1, fraction) &&
        all(isfinite, liquid) && all(>=(0), liquid) &&
        all(isfinite, ice) && all(>=(0), ice) &&
        all(isfinite, top) && all(>=(0), top) &&
        isapprox(sum(top), 1; atol=1e-12, rtol=0) ||
        error("finished cloud profile is invalid")
    return (; fraction, liquid, ice, top)
end

function _save_cloud_profile(capture::CloudProfileCapture, path::AbstractString)
    ispath(path) && error("cloud profile output already exists")
    data = _cloud_finish(capture)
    NCDatasets.NCDataset(path, "c") do ds
        NCDatasets.defDim(ds, "layer", length(capture.sigma))
        NCDatasets.defDim(ds, "cloud_top_bin", length(data.top))
        ds.attrib["schema"] = "readyesm_monthly_cloud_profile_v1"
        ds.attrib["start_time_UTC"] = string(capture.start)
        ds.attrib["stop_time_UTC"] = string(capture.stop)
        ds.attrib["state_sampling"] = "Hourly instants, start inclusive and stop exclusive"
        ds.attrib["state_samples"] = capture.samples
        ds.attrib["solver_order"] = "layer by point, bottom to top"
        ds.attrib["output_order"] = "point by layer, top to bottom before spatial averaging"
        ds.attrib["cloud_scheme"] = capture.scheme
        ds.attrib["cloud_top_proxy"] = "Topmost positive RRTMGP cloud fraction; column cover is maximum layer fraction"
        ds.attrib["cloud_top_clear_bin"] = length(data.top)
        for (name, values, dimension, units) in (
            ("sigma_full", capture.sigma, "layer", "1"),
            ("sigma_thickness", capture.dsigma, "layer", "1"),
            ("cloud_top_bin", Int32.(1:length(data.top)), "cloud_top_bin", "1"),
            ("cloud_fraction_profile", data.fraction, "layer", "1"),
            ("cloud_liquid_water_path_profile", data.liquid, "layer", "g m-2"),
            ("cloud_ice_water_path_profile", data.ice, "layer", "g m-2"),
            ("cloud_top_proxy_area_fraction", data.top, "cloud_top_bin", "1"))
            variable = NCDatasets.defVar(ds, name, eltype(values), (dimension,);
                attrib=Dict("units" => units))
            variable[:] = values
        end
    end
    return Dict("state_samples" => capture.samples, "file" => basename(path))
end
