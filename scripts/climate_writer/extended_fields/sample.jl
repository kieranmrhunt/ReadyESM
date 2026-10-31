# Included by the experimental runner after WriterIO and WriterSampling.
include("fields.jl")
using .LongRunFields
register_metadata!(WriterIO.WriterGeometry.SPECS)

function additional_samples(R, earth, ocean_weights, atmosphere_weights; ocean_latitude)
    O, RR = R.Oceananigans, R.RRTMGP
    atmosphere, ice = earth.atmosphere, earth.sea_ice.model
    p = atmosphere.variables.parameterizations
    radiation = atmosphere.model.longwave_radiation
    solver = radiation.solver
    npoints = length(atmosphere_weights)
    radiative = if radiation.call_counter > 0
        swd, swu = Array(RR.sw_flux_dn(solver)), Array(RR.sw_flux_up(solver))
        lwd, lwu = Array(RR.lw_flux_dn(solver)), Array(RR.lw_flux_up(solver))
        isnothing(solver.clear_flux_lw) && error("extended all-sky writer requires clear-sky fluxes")
        cswu, clwu = Array(RR.clear_sw_flux_up(solver)), Array(RR.clear_lw_flux_up(solver))
        radiation_fields(swd, swu, lwd, lwu, cswu, clwu)
    else
        earth.clock.iteration == 0 || error("radiation has no solve after a completed step")
        # A built model may not yet have solved radiation. Do not inspect raw
        # solver buffers or substitute zero as a physical launch-state flux.
        Dict(name => fill(NaN, npoints) for (name, (kind, unit, _)) in SPECS
            if kind == :atmosphere && unit == "W m-2" && name != "surface_sensible_heat_upward")
    end
    monthly = Dict{String,Any}(radiative)
    h = Array(O.interior(ice.ice_thickness))[:, :, 1]
    concentration = Array(O.interior(ice.ice_concentration))[:, :, 1]
    merge!(monthly, ice_fields(h, concentration))
    rain, snow = vec(Array(p.rain_rate.data)), vec(Array(p.snow_rate.data))
    monthly["snow_rate"] = Float64.(snow)
    monthly["total_precipitation_rate"] = Float64.(rain) .+ Float64.(snow)
    monthly["surface_sensible_heat_upward"] = Float64.(vec(Array(p.sensible_heat_flux.data)))
    monthly["surface_water_vapor_upward"] = Float64.(vec(Array(p.surface_humidity_flux.data)))
    mean(values, weights) = WriterSampling.weighted_mean(values, weights)
    # Launch radiation can be unavailable. Initial values are written separately
    # and never enter period accumulators. All later samples must be finite.
    toa = radiation.call_counter == 0 ? NaN : mean(monthly["toa_net_downward"], atmosphere_weights)
    daily = Dict(
        "global_toa_net_downward" => [toa],
        "global_total_precipitation_rate" => [mean(monthly["total_precipitation_rate"], atmosphere_weights)],
        "global_ice_area" => [mean(concentration, ocean_weights) * sum(ocean_weights)],
        "global_ice_volume" => [mean(monthly["ice_volume_per_area"], ocean_weights) * sum(ocean_weights)],
    )
    size(ocean_latitude) == size(ocean_weights) || throw(DimensionMismatch("ocean latitude and weights"))
    all(isfinite, ocean_latitude) || error("nonfinite ocean latitude")
    for (hemisphere, mask) in (("northern", ocean_latitude .>= 0), ("southern", ocean_latitude .< 0))
        weights = ocean_weights .* mask
        area = sum(weights)
        daily[hemisphere * "_ice_area"] = [iszero(area) ? 0. : mean(concentration, weights) * area]
        daily[hemisphere * "_ice_volume"] = [iszero(area) ? 0. : mean(monthly["ice_volume_per_area"], weights) * area]
    end
    return monthly, daily
end
