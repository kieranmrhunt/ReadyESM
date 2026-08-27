function _host_array(field)
    return Array(field.data)
end
_host_array(field::AbstractArray) = Array(field)

_netcdf_attribute_value(value) = value isa Bool ? Int8(value) : value

function _land_surface_field(field, land_fraction)
    host_field = _host_array(field)
    values = fill(NaN, length(land_fraction))
    for ij in eachindex(land_fraction)
        if land_fraction[ij] > 0
            values[ij] = Float64(_land_surface_temperature(host_field, ij))
        end
    end
    return values
end

function _column_cloud_fraction(model, npoints)
    radiation = model.longwave_radiation
    if hasproperty(radiation, :cloud_scheme) &&
       radiation.cloud_scheme in _ALL_SKY_CLOUD_SCHEMES
        layer_fraction = Array(RRTMGP.cloud_fraction(radiation.solver))
        return vec(Float64.(maximum(layer_fraction; dims = 1)))
    end
    return zeros(npoints)
end

"""Return RRTMGP cloud layers in SpeedyWeather point/top-to-bottom ordering."""
function _atmosphere_cloud_layer_diagnostics(model)
    npoints = model.geometry.npoints
    nlayers = model.geometry.nlayers
    radiation = model.longwave_radiation
    if !(radiation isa RRTMGPRadiation) ||
       !(radiation.cloud_scheme in _ALL_SKY_CLOUD_SCHEMES)
        empty_cloud = zeros(Float64, npoints, nlayers)
        return (
            fraction = empty_cloud,
            liquid_water_path = copy(empty_cloud),
            ice_water_path = copy(empty_cloud),
            column_fraction = zeros(Float64, npoints),
        )
    end

    solver_fraction = Float64.(Array(RRTMGP.cloud_fraction(radiation.solver)))
    solver_liquid = Float64.(Array(
        RRTMGP.cloud_liquid_water_path(radiation.solver),
    ))
    solver_ice = Float64.(Array(
        RRTMGP.cloud_ice_water_path(radiation.solver),
    ))
    expected_shape = (nlayers, npoints)
    size(solver_fraction) == expected_shape || throw(DimensionMismatch(
        "RRTMGP cloud-fraction shape $(size(solver_fraction)) does not match " *
        "expected layer/point shape $expected_shape",
    ))
    size(solver_liquid) == expected_shape && size(solver_ice) == expected_shape ||
        throw(DimensionMismatch("RRTMGP cloud condensate shapes differ"))

    # RRTMGP stores bottom-to-top layer by column; durable atmospheric output
    # follows SpeedyWeather's point by top-to-bottom layer convention.
    to_speedy_order(values) = reverse(permutedims(values, (2, 1)); dims = 2)
    fraction = to_speedy_order(solver_fraction)
    return (
        fraction,
        liquid_water_path = to_speedy_order(solver_liquid),
        ice_water_path = to_speedy_order(solver_ice),
        column_fraction = vec(maximum(fraction; dims = 2)),
    )
end

"""Return the surface temperature retained by the active radiation scheme."""
function _atmosphere_radiative_surface_temperature(vars, model)
    radiation = model.longwave_radiation
    if radiation isa RRTMGPRadiation
        values = Float64.(Array(RRTMGP.surface_temperature(radiation.solver)))
        length(values) == model.geometry.npoints || throw(DimensionMismatch(
            "RRTMGP surface-temperature length $(length(values)) does not " *
            "match $(model.geometry.npoints) atmospheric columns",
        ))
        return values
    end

    # Simplified SpeedyWeather radiation has no retained solver surface state.
    # Diagnose the same fourth-power mixture used by ReadyESM's RRTMGP adapter
    # so the output still has a physical, explicitly documented value.
    return Float64[
        _surface_temperature(ij, vars, model) for ij in 1:model.geometry.npoints
    ]
end

"""Return the surface albedo retained by the active shortwave solver."""
function _atmosphere_radiative_surface_albedo(vars, model)
    radiation = model.longwave_radiation
    if radiation isa RRTMGPRadiation
        direct = Float64.(Array(
            RRTMGP.direct_sw_surface_albedo(radiation.solver),
        ))
        diffuse = Float64.(Array(
            RRTMGP.diffuse_sw_surface_albedo(radiation.solver),
        ))
        size(direct) == size(diffuse) && ndims(direct) == 2 &&
            size(direct, 1) > 0 &&
            size(direct, 2) == model.geometry.npoints ||
            throw(DimensionMismatch(
                "RRTMGP direct/diffuse surface-albedo shapes " *
                "$(size(direct))/$(size(diffuse)) do not match " *
                "$(model.geometry.npoints) atmospheric columns",
            ))
        maximum(abs, direct .- diffuse) <= 8eps(Float32) || error(
            "RRTMGP direct and diffuse surface albedos differ",
        )
        albedo = vec(direct[1, :])
        maximum(abs, direct .- reshape(albedo, 1, :)) <= 8eps(Float32) ||
            error("RRTMGP surface albedo is not band invariant")
        all((0 .<= albedo) .& (albedo .<= 1)) || error(
            "RRTMGP surface albedo is outside [0, 1]",
        )
        return albedo
    end

    values = Float64.(Array(vars.parameterizations.albedo.data))
    length(values) == model.geometry.npoints || throw(DimensionMismatch(
        "SpeedyWeather surface-albedo length $(length(values)) does not " *
        "match $(model.geometry.npoints) atmospheric columns",
    ))
    albedo = vec(values)
    all(isfinite, albedo) && all((0 .<= albedo) .& (albedo .<= 1)) || error(
        "SpeedyWeather surface albedo is non-finite or outside [0, 1]",
    )
    return albedo
end

"""Machine-readable provenance for the radiative forcing actually in use."""
function _forcing_provenance(radiation, config::ExperimentConfig)
    rrtmgp_active = radiation isa RRTMGPRadiation
    aerosol_active = rrtmgp_active && radiation.target_aod_550nm > 0
    return (
        requested_co2_ppm = config.forcing.co2_ppm,
        applied_co2_ppm = rrtmgp_active ? Float64(radiation.co2_ppm) : NaN,
        co2_forcing_active = rrtmgp_active,
        co2_forcing_mode = rrtmgp_active ?
            "prescribed_well_mixed_volume_mixing_ratio" :
            "not_applied_by_speedy_simplified_radiation",
        requested_aerosol_optical_depth_550nm =
            config.forcing.aerosol_optical_depth_550nm,
        applied_aerosol_optical_depth_550nm = rrtmgp_active ?
            Float64(radiation.target_aod_550nm) : NaN,
        aerosol_species = config.forcing.aerosol_species,
        aerosol_forcing_active = aerosol_active,
        aerosol_forcing_mode = aerosol_active ?
            "prescribed_column_aod_550nm_recalibrated_each_radiation_solve" :
            "disabled",
        aerosol_vertical_profile = aerosol_active ?
            "fixed_log_pressure_gaussian_center_150hPa_width_0.8" :
            "none",
        forcing_time_dependence = "constant_for_entire_integration",
        interactive_carbon_cycle = false,
        interactive_aerosol_chemistry = false,
        aerosol_emissions = false,
        radiation_every_n_steps = config.forcing.radiation_every_n_steps,
        cloud_humidity_search_min_sigma =
            config.forcing.cloud_humidity_search_min_sigma,
        applied_cloud_humidity_search_min_sigma =
            rrtmgp_active &&
                config.forcing.cloud_scheme in _ALL_SKY_CLOUD_SCHEMES ?
            Float64(radiation.cloud_humidity_search_min_sigma) : NaN,
        cloud_humidity_search_scheme =
            config.forcing.cloud_scheme in _ALL_SKY_CLOUD_SCHEMES ?
            _diagnostic_cloud_humidity_search_scheme(
                config.forcing.cloud_humidity_search_min_sigma,
            ) : "disabled",
        cloud_condensate_retention_fraction =
            config.forcing.cloud_condensate_retention_fraction,
        cloud_condensate_residence_time_hours =
            config.forcing.cloud_condensate_residence_time_hours,
    )
end

"""Extract compact, architecture-independent diagnostics from a completed run."""
function collect_diagnostics(simulation, config::ExperimentConfig)
    vars = simulation.variables
    model = simulation.model
    grid = vars.prognostic.ocean.sea_surface_temperature.grid
    longitude, latitude = RG.get_londlatds(grid)
    land_fraction = Float64.(_host_array(model.land_sea_mask.mask))

    surface_temperature = model.callbacks[:global_surface_temperature].temperature
    radiation_budget = get(model.callbacks, :global_radiation_budget, nothing)
    timestep_days = Float64(model.time_stepping.Δt_sec) / 86_400
    time_days = collect((0:(length(surface_temperature) - 1)) .* timestep_days)

    return (
        longitude = Float64.(longitude),
        latitude = Float64.(latitude),
        time_days,
        global_surface_temperature = Float64.(surface_temperature),
        toa_incoming_shortwave = isnothing(radiation_budget) ? Float64[] :
            Float64.(radiation_budget.incoming_shortwave),
        toa_outgoing_shortwave = isnothing(radiation_budget) ? Float64[] :
            Float64.(radiation_budget.outgoing_shortwave),
        toa_outgoing_longwave = isnothing(radiation_budget) ? Float64[] :
            Float64.(radiation_budget.outgoing_longwave),
        toa_clear_outgoing_shortwave = isnothing(radiation_budget) ? Float64[] :
            Float64.(radiation_budget.clear_outgoing_shortwave),
        toa_clear_outgoing_longwave = isnothing(radiation_budget) ? Float64[] :
            Float64.(radiation_budget.clear_outgoing_longwave),
        toa_net_downward = isnothing(radiation_budget) ? Float64[] :
            Float64.(radiation_budget.net_downward),
        toa_clear_net_downward = isnothing(radiation_budget) ? Float64[] :
            Float64.(radiation_budget.clear_net_downward),
        surface_air_temperature = Float64.(
            _host_array(vars.parameterizations.surface_air_temperature),
        ),
        sea_surface_temperature = Float64.(
            _host_array(vars.prognostic.ocean.sea_surface_temperature),
        ),
        sea_ice_concentration = Float64.(
            _host_array(vars.prognostic.ocean.sea_ice_concentration),
        ),
        land_fraction,
        land_surface_temperature = _land_surface_field(
            vars.prognostic.land.soil_temperature,
            land_fraction,
        ),
        land_surface_moisture = _land_surface_field(
            vars.prognostic.land.soil_moisture,
            land_fraction,
        ),
        column_cloud_fraction = _column_cloud_fraction(model, length(longitude)),
        outgoing_longwave = Float64.(_host_array(vars.parameterizations.outgoing_longwave)),
        outgoing_shortwave = Float64.(_host_array(vars.parameterizations.outgoing_shortwave)),
        surface_longwave_down = Float64.(
            _host_array(vars.parameterizations.surface_longwave_down),
        ),
        surface_shortwave_down = Float64.(
            _host_array(vars.parameterizations.surface_shortwave_down),
        ),
        metadata = (
            experiment = config.name,
            truncation = config.truncation,
            nlayers = config.nlayers,
            atmosphere_hyperdiffusion_hours = config.atmosphere_hyperdiffusion_hours,
            atmosphere_divergence_hyperdiffusion_hours =
                config.atmosphere_divergence_hyperdiffusion_hours,
            atmosphere_vertical_diffusion =
                _atmosphere_vertical_diffusion_provenance(
                    model.vertical_diffusion,
                ),
            mixed_layer_depth_m = config.mixed_layer_depth_m,
            atmosphere_initial_conditions = String(config.atmosphere_initial_conditions),
            era5_month = config.era5_month,
            _atmosphere_initial_condition_provenance(
                model.initial_conditions,
                config,
            )...,
            land_model = String(config.land_model),
            _land_initial_condition_provenance(config)...,
            terrarium_soil_layers = config.terrarium_soil_layers,
            cloud_scheme = String(config.forcing.cloud_scheme),
            cloud_liquid_water_path_gm2 = config.forcing.cloud_liquid_water_path_gm2,
            cloud_ice_water_path_gm2 = config.forcing.cloud_ice_water_path_gm2,
            active_radiation = String(config.forcing.radiation),
            cloud_radiation_active =
                config.forcing.radiation == :rrtmgp_all_sky &&
                config.forcing.cloud_scheme in _ALL_SKY_CLOUD_SCHEMES,
            atmosphere_cloud_condensate_state =
                config.forcing.cloud_scheme == :prognostic_condensate ?
                "column_local_prognostic_liquid_ice_path_conservative_precipitation_delay_v1" :
                "not_prognostic_diagnostic_cloud_optics_only",
            _prognostic_cloud_condensate_summary(
                model.longwave_radiation,
                _global_point_weights(model.spectral_grid),
            )...,
            _forcing_provenance(model.longwave_radiation, config)...,
        ),
    )
end

function _write_netcdf(path::AbstractString, diagnostics)
    NCDataset(path, "c") do dataset
        defDim(dataset, "point", length(diagnostics.longitude))
        defDim(dataset, "time", length(diagnostics.time_days))

        longitude = defVar(dataset, "longitude", Float64, ("point",))
        latitude = defVar(dataset, "latitude", Float64, ("point",))
        time = defVar(dataset, "time", Float64, ("time",))
        longitude.attrib["units"] = "degrees_east"
        latitude.attrib["units"] = "degrees_north"
        time.attrib["units"] = "days"
        time.attrib["long_name"] = "elapsed model time"
        longitude[:] = diagnostics.longitude
        latitude[:] = diagnostics.latitude
        time[:] = diagnostics.time_days

        variables = (
            ("global_surface_temperature", diagnostics.global_surface_temperature, ("time",), "K"),
            ("surface_air_temperature", diagnostics.surface_air_temperature, ("point",), "K"),
            ("sea_surface_temperature", diagnostics.sea_surface_temperature, ("point",), "K"),
            ("sea_ice_concentration", diagnostics.sea_ice_concentration, ("point",), "1"),
            ("land_fraction", diagnostics.land_fraction, ("point",), "1"),
            ("land_surface_temperature", diagnostics.land_surface_temperature, ("point",), "K"),
            ("land_surface_moisture", diagnostics.land_surface_moisture, ("point",), "1"),
            ("column_cloud_fraction", diagnostics.column_cloud_fraction, ("point",), "1"),
            ("outgoing_longwave", diagnostics.outgoing_longwave, ("point",), "W/m^2"),
            ("outgoing_shortwave", diagnostics.outgoing_shortwave, ("point",), "W/m^2"),
            ("surface_longwave_down", diagnostics.surface_longwave_down, ("point",), "W/m^2"),
            ("surface_shortwave_down", diagnostics.surface_shortwave_down, ("point",), "W/m^2"),
        )
        for (name, values, dimensions, units) in variables
            variable = defVar(dataset, name, Float64, dimensions)
            variable.attrib["units"] = units
            variable[:] = values
        end

        if !isempty(diagnostics.toa_net_downward)
            budget_variables = (
                ("toa_incoming_shortwave", diagnostics.toa_incoming_shortwave),
                ("toa_outgoing_shortwave", diagnostics.toa_outgoing_shortwave),
                ("toa_outgoing_longwave", diagnostics.toa_outgoing_longwave),
                ("toa_clear_outgoing_shortwave", diagnostics.toa_clear_outgoing_shortwave),
                ("toa_clear_outgoing_longwave", diagnostics.toa_clear_outgoing_longwave),
                ("toa_net_downward", diagnostics.toa_net_downward),
                ("toa_clear_net_downward", diagnostics.toa_clear_net_downward),
            )
            for (name, values) in budget_variables
                variable = defVar(dataset, name, Float64, ("time",))
                variable.attrib["units"] = "W/m^2"
                variable[:] = values
            end
        end

        for (key, value) in pairs(diagnostics.metadata)
            # NetCDF attributes have no native boolean type.
            dataset.attrib[String(key)] = _netcdf_attribute_value(value)
        end
        dataset.attrib["created_at_utc"] = string(now(UTC))
        dataset.attrib["source"] = "ReadyESM coupled experiment"
    end
    return path
end

function _map_panel!(figure, position, title, diagnostics, values; colorrange = nothing)
    axis = Axis(
        _plot_layout(figure, position...);
        title,
        xlabel = "longitude (°E)",
        ylabel = "latitude (°N)",
        limits = ((0, 360), (-90, 90)),
    )
    plot = if isnothing(colorrange)
        scatter!(
            axis,
            diagnostics.longitude,
            diagnostics.latitude;
            color = values,
            markersize = 5,
        )
    else
        scatter!(
            axis,
            diagnostics.longitude,
            diagnostics.latitude;
            color = values,
            markersize = 5,
            colorrange,
        )
    end
    Colorbar(_plot_layout(figure, position[1], position[2] + 1), plot)
    return axis
end

function _write_figure(path::AbstractString, diagnostics)
    has_land_state = hasproperty(diagnostics, :land_surface_temperature)
    has_cloud_state = hasproperty(diagnostics, :column_cloud_fraction)
    has_budget = hasproperty(diagnostics, :toa_net_downward) && !isempty(diagnostics.toa_net_downward)
    figure = Figure(size = has_budget ? (1200, 1550) :
        (has_land_state || has_cloud_state ? (1200, 1300) : (1200, 800)))
    time_axis = Axis(
        _plot_layout(figure, 1, 1:4);
        title = "Global mean lowest-level atmospheric temperature",
        xlabel = "model time (days)",
        ylabel = "temperature (K)",
    )
    lines!(
        time_axis,
        diagnostics.time_days,
        diagnostics.global_surface_temperature;
        linewidth = 2,
    )

    map_start_row = 2
    if has_budget
        budget_axis = Axis(
            _plot_layout(figure, 2, 1:4);
            title = "Global top-of-atmosphere radiation budget",
            xlabel = "model time (days)",
            ylabel = "flux (W m⁻²)",
        )
        lines!(budget_axis, diagnostics.time_days, diagnostics.toa_incoming_shortwave;
               label = "incoming SW", linewidth = 2)
        lines!(budget_axis, diagnostics.time_days, diagnostics.toa_outgoing_shortwave;
               label = "reflected SW", linewidth = 2)
        lines!(budget_axis, diagnostics.time_days, diagnostics.toa_outgoing_longwave;
               label = "outgoing LW", linewidth = 2)
        lines!(budget_axis, diagnostics.time_days, diagnostics.toa_net_downward;
               label = "net downward", linewidth = 2)
        lines!(budget_axis, diagnostics.time_days, diagnostics.toa_clear_net_downward;
               label = "clear-sky net downward", linewidth = 2, linestyle = :dash)
        axislegend(budget_axis; position = :rb, orientation = :horizontal)
        map_start_row = 3
    end

    _map_panel!(
        figure,
        (map_start_row, 1),
        "Final surface air temperature (K)",
        diagnostics,
        diagnostics.surface_air_temperature,
    )
    _map_panel!(
        figure,
        (map_start_row, 3),
        "Final sea-surface temperature (K)",
        diagnostics,
        diagnostics.sea_surface_temperature,
    )
    _map_panel!(
        figure,
        (map_start_row + 1, 1),
        "Final sea-ice concentration",
        diagnostics,
        diagnostics.sea_ice_concentration;
        colorrange = (0, 1),
    )
    _map_panel!(
        figure,
        (map_start_row + 1, 3),
        "Land fraction",
        diagnostics,
        diagnostics.land_fraction;
        colorrange = (0, 1),
    )
    if has_land_state
        _map_panel!(
            figure,
            (map_start_row + 2, 1),
            "Final land-surface temperature (K)",
            diagnostics,
            diagnostics.land_surface_temperature,
        )
        _map_panel!(
            figure,
            (map_start_row + 2, 3),
            "Final land-surface soil moisture",
            diagnostics,
            diagnostics.land_surface_moisture;
            colorrange = (0, 1),
        )
    end
    if has_cloud_state
        _map_panel!(
            figure,
            (map_start_row + 3, 1),
            "Final column cloud fraction",
            diagnostics,
            diagnostics.column_cloud_fraction;
            colorrange = (0, 1),
        )
        if hasproperty(diagnostics, :outgoing_longwave)
            _map_panel!(
                figure,
                (map_start_row + 3, 3),
                "Final outgoing longwave (W m⁻²)",
                diagnostics,
                diagnostics.outgoing_longwave,
            )
        end
    end
    save(path, figure; px_per_unit = 1.5)
    return path
end

"""Render the summary PNG again from an existing ReadyESM diagnostic NetCDF file."""
function render_diagnostics(netcdf_path::AbstractString, figure_path::AbstractString)
    diagnostics = NCDataset(netcdf_path, "r") do dataset
        core = (
            longitude = Array(dataset["longitude"][:]),
            latitude = Array(dataset["latitude"][:]),
            time_days = Array(dataset["time"][:]),
            global_surface_temperature = Array(dataset["global_surface_temperature"][:]),
            surface_air_temperature = Array(dataset["surface_air_temperature"][:]),
            sea_surface_temperature = Array(dataset["sea_surface_temperature"][:]),
            sea_ice_concentration = Array(dataset["sea_ice_concentration"][:]),
            land_fraction = Array(dataset["land_fraction"][:]),
        )
        if haskey(dataset, "land_surface_temperature")
            core = merge(
                core,
                (
                    land_surface_temperature = Array(dataset["land_surface_temperature"][:]),
                    land_surface_moisture = Array(dataset["land_surface_moisture"][:]),
                ),
            )
        end
        if haskey(dataset, "column_cloud_fraction")
            core = merge(
                core,
                (
                    column_cloud_fraction = Array(dataset["column_cloud_fraction"][:]),
                    outgoing_longwave = Array(dataset["outgoing_longwave"][:]),
                ),
            )
        end
        if haskey(dataset, "toa_net_downward")
            core = merge(
                core,
                (
                    toa_incoming_shortwave = Array(dataset["toa_incoming_shortwave"][:]),
                    toa_outgoing_shortwave = Array(dataset["toa_outgoing_shortwave"][:]),
                    toa_outgoing_longwave = Array(dataset["toa_outgoing_longwave"][:]),
                    toa_net_downward = Array(dataset["toa_net_downward"][:]),
                    toa_clear_net_downward = haskey(dataset, "toa_clear_net_downward") ?
                        Array(dataset["toa_clear_net_downward"][:]) :
                        Array(dataset["toa_net_downward"][:]),
                ),
            )
        end
        core
    end
    mkpath(dirname(figure_path))
    return _write_figure(figure_path, diagnostics)
end

"""Write a NetCDF diagnostic bundle and a PNG summary figure."""
function save_diagnostics(diagnostics, config::ExperimentConfig)
    mkpath(config.output_dir)
    netcdf_path = joinpath(config.output_dir, "diagnostics.nc")
    figure_path = joinpath(config.output_dir, "summary.png")
    _write_netcdf(netcdf_path, diagnostics)
    _write_figure(figure_path, diagnostics)
    return (; netcdf_path, figure_path)
end
"""Time-weighted means in complete, end-aligned fixed-length annual blocks."""
function _annual_block_means(
    time::AbstractVector,
    values::AbstractVector,
    start_day,
    stop_day;
    year_days = 365.0,
)
    length(time) == length(values) || throw(DimensionMismatch(
        "annual-block time and value lengths differ",
    ))
    all(isfinite, time) && all(isfinite, values) || error(
        "annual-block input contains non-finite values",
    )
    all(>(0), diff(time)) || error(
        "annual-block time coordinate is not strictly increasing",
    )
    year_days > 0 || error("annual-block year length must be positive")
    stop_day > start_day || error("annual-block interval must increase")
    block_count = floor(Int, (stop_day - start_day) / year_days)
    block_count > 0 || error("annual-block interval contains no complete year")
    aligned_start = stop_day - block_count * year_days
    centers = Vector{Float64}(undef, block_count)
    means = Vector{Float64}(undef, block_count)
    tolerance = 64 * eps(Float64) * max(abs(start_day), abs(stop_day), 1.0)
    for block in 1:block_count
        left = aligned_start + (block - 1) * year_days
        right = left + year_days
        selected = (time .>= left - tolerance) .& (time .<= right + tolerance)
        block_time = Float64.(time[selected])
        block_values = Float64.(values[selected])
        length(block_time) >= 2 || error(
            "annual block $block has fewer than two samples",
        )
        abs(first(block_time) - left) <= tolerance || error(
            "annual block $block begins at $(first(block_time)), expected $left",
        )
        abs(last(block_time) - right) <= tolerance || error(
            "annual block $block ends at $(last(block_time)), expected $right",
        )
        intervals = diff(block_time)
        means[block] = sum(
            intervals .* (block_values[1:end-1] .+ block_values[2:end]) ./ 2,
        ) / sum(intervals)
        centers[block] = (left + right) / 2
    end
    return (; centers, means, aligned_start, year_days)
end

function _annual_block_means(
    time::AbstractVector,
    values::AbstractMatrix,
    start_day,
    stop_day;
    year_days = 365.0,
)
    size(values, 2) == length(time) || throw(DimensionMismatch(
        "annual-block profile and time lengths differ",
    ))
    size(values, 1) > 0 || error("annual-block profile has no levels")
    all(isfinite, time) && all(isfinite, values) || error(
        "annual-block profile input contains non-finite values",
    )
    all(>(0), diff(time)) || error(
        "annual-block time coordinate is not strictly increasing",
    )
    year_days > 0 || error("annual-block year length must be positive")
    stop_day > start_day || error("annual-block interval must increase")
    block_count = floor(Int, (stop_day - start_day) / year_days)
    block_count > 0 || error("annual-block interval contains no complete year")
    aligned_start = stop_day - block_count * year_days
    centers = Vector{Float64}(undef, block_count)
    means = Matrix{Float64}(undef, size(values, 1), block_count)
    tolerance = 64 * eps(Float64) * max(abs(start_day), abs(stop_day), 1.0)
    for block in 1:block_count
        left = aligned_start + (block - 1) * year_days
        right = left + year_days
        selected = (time .>= left - tolerance) .& (time .<= right + tolerance)
        block_time = Float64.(time[selected])
        block_values = Float64.(values[:, selected])
        length(block_time) >= 2 || error(
            "annual profile block $block has fewer than two samples",
        )
        abs(first(block_time) - left) <= tolerance || error(
            "annual profile block $block begins at $(first(block_time)), expected $left",
        )
        abs(last(block_time) - right) <= tolerance || error(
            "annual profile block $block ends at $(last(block_time)), expected $right",
        )
        intervals = diff(block_time)
        weighted = (
            block_values[:, 1:end-1] .+ block_values[:, 2:end]
        ) .* reshape(intervals, 1, :)
        means[:, block] .= vec(sum(weighted; dims = 2)) ./ (2 * sum(intervals))
        centers[block] = (left + right) / 2
    end
    return (; centers, means, aligned_start, year_days)
end

"""End-aligned complete Gregorian-year boundaries in elapsed model days."""
function _calendar_year_boundaries(
    model_start::DateTime,
    start_day,
    stop_day,
)
    all(isfinite, (start_day, stop_day)) || error(
        "calendar-year interval contains a non-finite day",
    )
    isinteger(start_day) && isinteger(stop_day) || error(
        "calendar-year boundaries require whole elapsed days",
    )
    stop_day > start_day || error("calendar-year interval must increase")
    analysis_date = model_start + Day(round(Int, start_day))
    final_date = model_start + Day(round(Int, stop_day))
    dates = DateTime[final_date]
    cursor = final_date
    while true
        previous = cursor - Year(1)
        previous < analysis_date && break
        push!(dates, previous)
        cursor = previous
    end
    reverse!(dates)
    length(dates) >= 2 || error(
        "calendar-year interval contains no complete Gregorian year",
    )
    milliseconds_per_day = 86_400_000
    boundaries = Float64[
        Dates.value(date - model_start) / milliseconds_per_day for date in dates
    ]
    all(>(0), diff(boundaries)) || error(
        "calendar-year boundaries are not strictly increasing",
    )
    return (; boundaries, dates)
end

function _calendar_annual_block_means(
    time::AbstractVector,
    values::AbstractVector,
    model_start::DateTime,
    start_day,
    stop_day,
)
    calendar = _calendar_year_boundaries(model_start, start_day, stop_day)
    block_count = length(calendar.boundaries) - 1
    centers = Vector{Float64}(undef, block_count)
    means = Vector{Float64}(undef, block_count)
    for block in 1:block_count
        left = calendar.boundaries[block]
        right = calendar.boundaries[block + 1]
        annual = _annual_block_means(
            time,
            values,
            left,
            right;
            year_days = right - left,
        )
        centers[block] = only(annual.centers)
        means[block] = only(annual.means)
    end
    return (
        ; centers,
        means,
        aligned_start = first(calendar.boundaries),
        boundary_days = calendar.boundaries,
        boundary_dates = calendar.dates,
    )
end

function _calendar_annual_block_means(
    time::AbstractVector,
    values::AbstractMatrix,
    model_start::DateTime,
    start_day,
    stop_day,
)
    calendar = _calendar_year_boundaries(model_start, start_day, stop_day)
    block_count = length(calendar.boundaries) - 1
    centers = Vector{Float64}(undef, block_count)
    means = Matrix{Float64}(undef, size(values, 1), block_count)
    for block in 1:block_count
        left = calendar.boundaries[block]
        right = calendar.boundaries[block + 1]
        annual = _annual_block_means(
            time,
            values,
            left,
            right;
            year_days = right - left,
        )
        centers[block] = only(annual.centers)
        means[:, block] .= vec(annual.means)
    end
    return (
        ; centers,
        means,
        aligned_start = first(calendar.boundaries),
        boundary_days = calendar.boundaries,
        boundary_dates = calendar.dates,
    )
end

"""Linear trend over the last complete annual means, expressed per decade."""
function _late_annual_trend_per_decade(
    centers::AbstractVector,
    annual_means::AbstractVector;
    years = 10,
    year_days = 365.2425,
)
    length(centers) == length(annual_means) || throw(DimensionMismatch(
        "annual centers and means differ in length",
    ))
    years >= 2 || error("annual trend requires at least two years")
    length(centers) >= years || error(
        "annual trend requires $years blocks, found $(length(centers))",
    )
    selected = (length(centers) - years + 1):length(centers)
    time = Float64.(centers[selected])
    values = Float64.(annual_means[selected])
    anomaly = time .- mean(time)
    denominator = sum(abs2, anomaly)
    denominator > 0 || error("annual trend time variance is zero")
    slope_per_day = sum(anomaly .* (values .- mean(values))) / denominator
    return slope_per_day * 10 * year_days
end

function _late_annual_trend_per_decade(
    centers::AbstractVector,
    annual_means::AbstractMatrix;
    years = 10,
    year_days = 365.2425,
)
    size(annual_means, 2) == length(centers) || throw(DimensionMismatch(
        "annual profile centers and means differ in length",
    ))
    return [
        _late_annual_trend_per_decade(
            centers,
            view(annual_means, level, :);
            years,
            year_days,
        ) for level in axes(annual_means, 1)
    ]
end
