struct ERA5Atmosphere{V, A3, A2} <: SpeedyWeather.AbstractInitialConditions
    pressure_levels_hpa::V
    latitude::V
    longitude::V
    temperature::A3
    zonal_wind::A3
    meridional_wind::A3
    specific_humidity::A3
    surface_pressure::A2
    surface_temperature::A2
    total_column_water_vapour::A2
    month::Int
    source_time::DateTime
    source_time_index::Int
    source_mode::Symbol
    pressure_levels_path::String
    single_levels_path::String
    pressure_levels_sha256::String
    single_levels_sha256::String
end

function _era5_array(variable, indices...)
    values = Array(variable[indices...])
    return Float32.(coalesce.(values, NaN32))
end

function _era5_sha256(path)
    return open(path, "r") do io
        bytes2hex(SHA.sha256(io))
    end
end

function _era5_normalized_units(variable)
    units = haskey(variable.attrib, "units") ? String(variable.attrib["units"]) : ""
    return replace(lowercase(units), " " => "", "**" => "^", "_" => "")
end

function _era5_require_units(variable, name, accepted)
    units = _era5_normalized_units(variable)
    units in accepted || error(
        "ERA5 $name units '$units' are not one of $(join(accepted, ", "))",
    )
    return nothing
end

function _era5_time_coordinate(dataset)
    name = if haskey(dataset, "valid_time")
        "valid_time"
    elseif haskey(dataset, "time")
        "time"
    else
        error("ERA5 file has neither valid_time nor time coordinate")
    end
    times = DateTime.(Array(dataset[name][:]))
    isempty(times) && error("ERA5 time coordinate is empty")
    issorted(times) || error("ERA5 time coordinate must increase")
    allunique(times) || error("ERA5 time coordinate contains duplicate timestamps")
    return name, times
end

function _era5_selected_time(dataset, config)
    _, times = _era5_time_coordinate(dataset)
    if config.atmosphere_initial_conditions == :era5_instantaneous
        matches = findall(==(config.start_date), times)
        length(matches) == 1 || error(
            "ERA5 instantaneous file must contain the model start timestamp " *
            "$(config.start_date) exactly once; found $(length(matches)) matches",
        )
    else
        matches = findall(times) do time
            year(time) == year(config.start_date) && month(time) == config.era5_month
        end
        length(matches) == 1 || error(
            "ERA5 monthly file must contain exactly one entry for " *
            "$(year(config.start_date))-$(lpad(config.era5_month, 2, '0')); " *
            "found $(length(matches)) matches",
        )
    end
    index = only(matches)
    return index, times[index], times
end

function _era5_require_instantaneous(variable, name)
    haskey(variable.attrib, "GRIB_stepType") || error(
        "ERA5 instantaneous $name lacks GRIB_stepType provenance",
    )
    String(variable.attrib["GRIB_stepType"]) == "instant" || error(
        "ERA5 instantaneous $name is marked GRIB_stepType=" *
        "$(variable.attrib["GRIB_stepType"]), not instant",
    )
    return nothing
end

function _era5_validate_finite_range(field, name, lower, upper; missing_allowed = false)
    values = Float64.(vec(field))
    finite = values[isfinite.(values)]
    isempty(finite) && error("ERA5 $name contains no finite values")
    !missing_allowed && length(finite) != length(values) && error(
        "ERA5 $name contains non-finite values",
    )
    extrema(finite)[1] >= lower || error("ERA5 $name is below $lower")
    extrema(finite)[2] <= upper || error("ERA5 $name exceeds $upper")
    return nothing
end

function ERA5Atmosphere(config::ExperimentConfig)
    month = config.era5_month
    source_mode = config.atmosphere_initial_conditions == :era5_instantaneous ?
        :instantaneous : :monthly_mean
    pressure_levels_hpa = Float32[]
    latitude = Float32[]
    longitude = Float32[]
    temperature = Array{Float32}(undef, 0, 0, 0)
    zonal_wind = similar(temperature)
    meridional_wind = similar(temperature)
    specific_humidity = similar(temperature)
    pressure_time_index = 0
    source_time = DateTime(0)
    pressure_times = DateTime[]
    NCDataset(config.era5_pressure_levels_path, "r") do dataset
        pressure_time_index, source_time, pressure_times =
            _era5_selected_time(dataset, config)
        _era5_require_units(dataset["pressure_level"], "pressure coordinate", ("hpa",))
        _era5_require_units(dataset["latitude"], "latitude", ("degreesnorth",))
        _era5_require_units(dataset["longitude"], "longitude", ("degreeseast",))
        _era5_require_units(dataset["t"], "temperature", ("k",))
        _era5_require_units(dataset["u"], "zonal wind", ("ms^-1", "ms-1"))
        _era5_require_units(dataset["v"], "meridional wind", ("ms^-1", "ms-1"))
        _era5_require_units(
            dataset["q"],
            "specific humidity",
            ("kgkg^-1", "kgkg-1", "1"),
        )
        if source_mode == :instantaneous
            for (name, description) in (
                ("t", "temperature"),
                ("u", "zonal wind"),
                ("v", "meridional wind"),
                ("q", "specific humidity"),
            )
                _era5_require_instantaneous(dataset[name], description)
            end
        end
        pressure_levels_hpa = _era5_array(dataset["pressure_level"], :)
        latitude = _era5_array(dataset["latitude"], :)
        longitude = _era5_array(dataset["longitude"], :)
        # NCDatasets exposes the NetCDF dimensions in Julia order as
        # longitude, latitude, pressure, time. Store pressure, latitude,
        # longitude for the interpolation routines below.
        temperature = permutedims(
            _era5_array(dataset["t"], :, :, :, pressure_time_index),
            (3, 2, 1),
        )
        zonal_wind = permutedims(
            _era5_array(dataset["u"], :, :, :, pressure_time_index),
            (3, 2, 1),
        )
        meridional_wind = permutedims(
            _era5_array(dataset["v"], :, :, :, pressure_time_index),
            (3, 2, 1),
        )
        specific_humidity = permutedims(
            _era5_array(dataset["q"], :, :, :, pressure_time_index),
            (3, 2, 1),
        )
    end

    surface_pressure = Array{Float32}(undef, 0, 0)
    surface_temperature = similar(surface_pressure)
    total_column_water_vapour = similar(surface_pressure)
    NCDataset(config.era5_single_levels_path, "r") do dataset
        single_time_index, single_source_time, single_times =
            _era5_selected_time(dataset, config)
        single_source_time == source_time || error(
            "ERA5 pressure- and single-level selected timestamps differ",
        )
        single_times == pressure_times || error(
            "ERA5 pressure- and single-level time coordinates differ",
        )
        _era5_require_units(dataset["latitude"], "latitude", ("degreesnorth",))
        _era5_require_units(dataset["longitude"], "longitude", ("degreeseast",))
        _era5_require_units(dataset["sp"], "surface pressure", ("pa",))
        _era5_require_units(dataset["t2m"], "2-m temperature", ("k",))
        _era5_require_units(
            dataset["tcwv"],
            "total-column water vapour",
            ("kgm^-2", "kgm-2"),
        )
        if source_mode == :instantaneous
            for (name, description) in (
                ("sp", "surface pressure"),
                ("t2m", "2-m temperature"),
                ("tcwv", "total-column water vapour"),
            )
                _era5_require_instantaneous(dataset[name], description)
            end
        end
        surface_pressure = permutedims(
            _era5_array(dataset["sp"], :, :, single_time_index),
            (2, 1),
        )
        surface_temperature = permutedims(
            _era5_array(dataset["t2m"], :, :, single_time_index),
            (2, 1),
        )
        total_column_water_vapour = permutedims(
            _era5_array(dataset["tcwv"], :, :, single_time_index),
            (2, 1),
        )
        source_latitude = _era5_array(dataset["latitude"], :)
        source_longitude = _era5_array(dataset["longitude"], :)
        source_latitude == latitude || error(
            "ERA5 pressure- and single-level latitude coordinates differ",
        )
        source_longitude == longitude || error(
            "ERA5 pressure- and single-level longitude coordinates differ",
        )
    end

    nlevel = length(pressure_levels_hpa)
    nlatitude = length(latitude)
    nlongitude = length(longitude)
    expected_3d = (nlevel, nlatitude, nlongitude)
    size(temperature) == expected_3d || error(
        "unexpected ERA5 pressure-level shape $(size(temperature)); expected $expected_3d",
    )
    all(size(field) == expected_3d for field in (
        zonal_wind,
        meridional_wind,
        specific_humidity,
    )) || error("ERA5 pressure-level variables do not share one shape")
    expected_2d = (nlatitude, nlongitude)
    all(size(field) == expected_2d for field in (
        surface_pressure,
        surface_temperature,
        total_column_water_vapour,
    )) || error("ERA5 single-level variables do not share the coordinate shape")
    all(diff(longitude) .> 0) || error("ERA5 longitude must increase")
    all(diff(latitude) .< 0) || error("ERA5 latitude must decrease")
    all(diff(pressure_levels_hpa) .< 0) || error(
        "ERA5 pressure levels must be stored from high to low pressure",
    )
    if source_mode == :instantaneous
        nlevel >= 12 || error(
            "production ERA5 initialization requires at least 12 pressure levels; " *
            "found $nlevel",
        )
        maximum(pressure_levels_hpa) >= 1000 || error(
            "production ERA5 initialization must reach at least 1000 hPa",
        )
        minimum(pressure_levels_hpa) <= 100 || error(
            "production ERA5 initialization must reach at least 100 hPa",
        )
        # Require dense sampling through the troposphere, where interpolation
        # controls most mass and water. ERA5's standard stratospheric levels
        # intentionally widen above 100 hPa (for example 50 -> 30 and
        # 20 -> 10 hPa); rejecting those levels forces the model top to be
        # extrapolated from the troposphere instead of initialized from data.
        troposphere = pressure_levels_hpa[pressure_levels_hpa .>= 100]
        maximum(troposphere[1:end-1] ./ troposphere[2:end]) <= 1.6 || error(
            "production ERA5 tropospheric pressure-level spacing contains a " *
            "ratio larger than 1.6",
        )
    end

    _era5_validate_finite_range(temperature, "temperature", 150, 350; missing_allowed = true)
    _era5_validate_finite_range(zonal_wind, "zonal wind", -200, 200; missing_allowed = true)
    _era5_validate_finite_range(
        meridional_wind,
        "meridional wind",
        -200,
        200;
        missing_allowed = true,
    )
    _era5_validate_finite_range(
        specific_humidity,
        "specific humidity",
        0,
        0.1;
        missing_allowed = true,
    )
    _era5_validate_finite_range(surface_pressure, "surface pressure", 40_000, 115_000)
    _era5_validate_finite_range(surface_temperature, "2-m temperature", 180, 340)
    _era5_validate_finite_range(
        total_column_water_vapour,
        "total-column water vapour",
        0,
        100,
    )

    return ERA5Atmosphere(
        pressure_levels_hpa,
        latitude,
        longitude,
        temperature,
        zonal_wind,
        meridional_wind,
        specific_humidity,
        surface_pressure,
        surface_temperature,
        total_column_water_vapour,
        month,
        source_time,
        pressure_time_index,
        source_mode,
        config.era5_pressure_levels_path,
        config.era5_single_levels_path,
        _era5_sha256(config.era5_pressure_levels_path),
        _era5_sha256(config.era5_single_levels_path),
    )
end

@inline function _era5_horizontal_indices(ic, longitude, latitude)
    nlongitude = length(ic.longitude)
    longitude_spacing = ic.longitude[2] - ic.longitude[1]
    wrapped_longitude = mod(longitude - ic.longitude[1], 360)
    longitude_position = wrapped_longitude / longitude_spacing
    longitude_floor = floor(Int, longitude_position)
    i0 = mod(longitude_floor, nlongitude) + 1
    i1 = mod(i0, nlongitude) + 1
    longitude_fraction = longitude_position - longitude_floor

    nlatitude = length(ic.latitude)
    latitude_spacing = ic.latitude[2] - ic.latitude[1]
    clipped_latitude = clamp(latitude, ic.latitude[end], ic.latitude[1])
    latitude_position =
        clamp((clipped_latitude - ic.latitude[1]) / latitude_spacing, 0, nlatitude - 1)
    latitude_floor = floor(Int, latitude_position)
    j0 = latitude_floor + 1
    j1 = min(j0 + 1, nlatitude)
    latitude_fraction = j0 == j1 ? zero(latitude_position) :
        latitude_position - latitude_floor

    return i0, i1, longitude_fraction, j0, j1, latitude_fraction
end

@inline function _era5_weighted_valid(values, weights)
    numerator = 0.0
    denominator = 0.0
    for index in eachindex(values)
        value = values[index]
        weight = weights[index]
        if isfinite(value) && weight > 0
            numerator += weight * value
            denominator += weight
        end
    end
    return denominator > 0 ? numerator / denominator : NaN
end

@inline function _era5_bilinear(
    field::AbstractMatrix,
    ic,
    longitude,
    latitude,
)
    i0, i1, fi, j0, j1, fj = _era5_horizontal_indices(ic, longitude, latitude)
    values = (
        field[j0, i0],
        field[j0, i1],
        field[j1, i0],
        field[j1, i1],
    )
    weights = (
        (1 - fi) * (1 - fj),
        fi * (1 - fj),
        (1 - fi) * fj,
        fi * fj,
    )
    return _era5_weighted_valid(values, weights)
end

@inline function _era5_bilinear_level(
    field::AbstractArray{<:Any, 3},
    level,
    ic::ERA5Atmosphere,
    longitude,
    latitude,
)
    i0, i1, fi, j0, j1, fj = _era5_horizontal_indices(ic, longitude, latitude)
    values = (
        field[level, j0, i0],
        field[level, j0, i1],
        field[level, j1, i0],
        field[level, j1, i1],
    )
    weights = (
        (1 - fi) * (1 - fj),
        fi * (1 - fj),
        (1 - fi) * fj,
        fi * fj,
    )
    return _era5_weighted_valid(values, weights)
end

function _era5_profile_value(
    field,
    ic::ERA5Atmosphere,
    longitude,
    latitude,
    pressure_hpa;
    surface_pressure_hpa = Inf,
    surface_value = NaN,
    logarithmic_value = false,
)
    pressures = Float64[]
    values = Float64[]
    if isfinite(surface_value) && isfinite(surface_pressure_hpa)
        push!(pressures, surface_pressure_hpa)
        push!(values, surface_value)
    end
    for level in eachindex(ic.pressure_levels_hpa)
        level_pressure = ic.pressure_levels_hpa[level]
        value = _era5_bilinear_level(field, level, ic, longitude, latitude)
        isfinite(value) || continue
        push!(pressures, level_pressure)
        push!(values, value)
    end
    isempty(values) && return NaN
    length(values) == 1 && return values[1]

    order = sortperm(pressures)
    pressures = pressures[order]
    values = values[order]
    target = log(max(Float64(pressure_hpa), 1.0))
    coordinate = log.(pressures)
    upper = if target <= coordinate[1]
        2
    elseif target >= coordinate[end]
        length(coordinate)
    else
        searchsortedfirst(coordinate, target)
    end
    lower = upper - 1
    fraction = (target - coordinate[lower]) /
               (coordinate[upper] - coordinate[lower])
    if logarithmic_value
        lower_value = log(max(values[lower], 1e-12))
        upper_value = log(max(values[upper], 1e-12))
        return exp(lower_value + fraction * (upper_value - lower_value))
    end
    return values[lower] + fraction * (values[upper] - values[lower])
end

@inline function _era5_surface_pressure_hpa(ic, longitude, latitude)
    return _era5_bilinear(ic.surface_pressure, ic, longitude, latitude) / 100
end

@inline function _era5_temperature(ic, longitude, latitude, pressure_hpa)
    surface_pressure_hpa = _era5_surface_pressure_hpa(ic, longitude, latitude)
    surface_temperature =
        _era5_bilinear(ic.surface_temperature, ic, longitude, latitude)
    return clamp(_era5_profile_value(
        ic.temperature,
        ic,
        longitude,
        latitude,
        pressure_hpa;
        surface_pressure_hpa,
        surface_value = surface_temperature,
    ), 180, 330)
end

@inline function _era5_wind(field, ic, longitude, latitude, pressure_hpa)
    surface_pressure_hpa = _era5_surface_pressure_hpa(ic, longitude, latitude)
    return clamp(_era5_profile_value(
        field,
        ic,
        longitude,
        latitude,
        pressure_hpa;
        surface_pressure_hpa,
    ), -150, 150)
end

@inline function _era5_raw_humidity(ic, longitude, latitude, pressure_hpa)
    surface_pressure_hpa = _era5_surface_pressure_hpa(ic, longitude, latitude)
    return max(0, _era5_profile_value(
        ic.specific_humidity,
        ic,
        longitude,
        latitude,
        pressure_hpa;
        surface_pressure_hpa,
        logarithmic_value = true,
    ))
end

function _era5_humidity_scale(
    ic,
    longitude,
    latitude,
    sigma_full,
    sigma_half,
    atmosphere,
    gravity,
)
    surface_pressure_pa =
        _era5_bilinear(ic.surface_pressure, ic, longitude, latitude)
    target_tcwv =
        max(0, _era5_bilinear(ic.total_column_water_vapour, ic, longitude, latitude))
    raw_humidity = similar(sigma_full, Float64)
    saturation_limit = similar(sigma_full, Float64)
    for k in eachindex(sigma_full)
        pressure_pa = sigma_full[k] * surface_pressure_pa
        raw_humidity[k] =
            _era5_raw_humidity(ic, longitude, latitude, pressure_pa / 100)
        temperature = _era5_temperature(
            ic,
            longitude,
            latitude,
            pressure_pa / 100,
        )
        saturation_limit[k] = 0.98 * SpeedyWeather.saturation_humidity(
            temperature,
            pressure_pa,
            atmosphere,
        )
    end

    column_water(scale) = surface_pressure_pa / gravity * sum(
        min(scale * raw_humidity[k], saturation_limit[k]) *
        (sigma_half[k + 1] - sigma_half[k]) for k in eachindex(sigma_full)
    )
    target_tcwv <= 0 && return 0.0
    upper = max(1.0, target_tcwv / max(column_water(1.0), eps(Float64)))
    while column_water(upper) < target_tcwv && upper < 65_536
        upper *= 2
    end
    lower = 0.0
    for _ in 1:24
        middle = (lower + upper) / 2
        if column_water(middle) < target_tcwv
            lower = middle
        else
            upper = middle
        end
    end
    return (lower + upper) / 2
end

function SpeedyWeather.initialize!(
    vars::SpeedyWeather.Variables,
    initial_conditions::ERA5Atmosphere,
    model::SpeedyWeather.PrimitiveWet,
)
    sigma_full = Float64.(Array(model.geometry.σ_levels_full))
    sigma_half = Float64.(Array(model.geometry.σ_levels_half))
    atmosphere = model.atmosphere
    gravity = Float64(model.planet.gravity)

    pressure_function = (longitude, latitude) -> log(
        100 * _era5_surface_pressure_hpa(initial_conditions, longitude, latitude),
    )
    temperature_function = (longitude, latitude, sigma) -> begin
        surface_pressure_hpa =
            _era5_surface_pressure_hpa(initial_conditions, longitude, latitude)
        _era5_temperature(
            initial_conditions,
            longitude,
            latitude,
            sigma * surface_pressure_hpa,
        )
    end
    zonal_wind_function = (longitude, latitude, sigma) -> begin
        surface_pressure_hpa =
            _era5_surface_pressure_hpa(initial_conditions, longitude, latitude)
        _era5_wind(
            initial_conditions.zonal_wind,
            initial_conditions,
            longitude,
            latitude,
            sigma * surface_pressure_hpa,
        )
    end
    meridional_wind_function = (longitude, latitude, sigma) -> begin
        surface_pressure_hpa =
            _era5_surface_pressure_hpa(initial_conditions, longitude, latitude)
        _era5_wind(
            initial_conditions.meridional_wind,
            initial_conditions,
            longitude,
            latitude,
            sigma * surface_pressure_hpa,
        )
    end
    humidity_function = (longitude, latitude, sigma) -> begin
        surface_pressure_pa =
            _era5_bilinear(
                initial_conditions.surface_pressure,
                initial_conditions,
                longitude,
                latitude,
            )
        pressure_pa = sigma * surface_pressure_pa
        scale = _era5_humidity_scale(
            initial_conditions,
            longitude,
            latitude,
            sigma_full,
            sigma_half,
            atmosphere,
            gravity,
        )
        raw_humidity = _era5_raw_humidity(
            initial_conditions,
            longitude,
            latitude,
            pressure_pa / 100,
        )
        temperature = _era5_temperature(
            initial_conditions,
            longitude,
            latitude,
            pressure_pa / 100,
        )
        saturation_limit = 0.98 * SpeedyWeather.saturation_humidity(
            temperature,
            pressure_pa,
            atmosphere,
        )
        min(scale * raw_humidity, saturation_limit)
    end

    # SpeedyWeather's three-dimensional function setter honours
    # `static_func=false` by evaluating dynamic closures on the CPU. Its
    # two-dimensional setter does not: it broadcasts the closure directly on
    # the field architecture. ERA5 owns ordinary host arrays and strings, so
    # that closure cannot be a CUDA kernel argument. Evaluate log surface
    # pressure on a host RingGrid explicitly, transfer only the numeric field,
    # and then perform the spectral transform on the model architecture.
    NF = model.spectral_grid.NF
    grid = model.spectral_grid.grid
    pressure_grid = zeros(NF, grid)
    pressure_grid_cpu = SpeedyWeather.on_architecture(
        SpeedyWeather.CPU(),
        pressure_grid,
    )
    SpeedyWeather.set!(pressure_grid_cpu, pressure_function)
    pressure_grid = SpeedyWeather.on_architecture(
        model.spectral_grid.architecture,
        pressure_grid_cpu,
    )
    pressure = SpeedyWeather.get_step(vars.prognostic.pressure, 1)
    SpeedyWeather.set!(
        pressure,
        pressure_grid,
        model.geometry,
        model.spectral_transform,
    )
    SpeedyWeather.set!(
        vars,
        model;
        temperature = temperature_function,
        lf = 1,
        static_func = false,
    )
    # SpeedyWeather's set_vordiv!(::LowerTriangularArray, ...) implementation
    # accidentally constructs its meridional spectral field from `u` again.
    # Build both grid fields through the standard setter, then apply the same
    # radius-aware curl/divergence that the public route intends. This stores
    # physical s⁻¹ vorticity/divergence; run! performs the model's later radius
    # scaling exactly once.
    nlayers = model.geometry.nlayers
    zonal_wind_grid = zeros(NF, grid, nlayers)
    meridional_wind_grid = zeros(NF, grid, nlayers)
    SpeedyWeather.set!(
        zonal_wind_grid,
        zonal_wind_function,
        model.geometry,
        model.spectral_transform;
        static_func = false,
    )
    SpeedyWeather.set!(
        meridional_wind_grid,
        meridional_wind_function,
        model.geometry,
        model.spectral_transform;
        static_func = false,
    )
    zonal_scaled = SpeedyWeather.RingGrids.scale_coslat⁻¹(zonal_wind_grid)
    meridional_scaled = SpeedyWeather.RingGrids.scale_coslat⁻¹(meridional_wind_grid)
    zonal_spectral = SpeedyWeather.SpeedyTransforms.transform(
        zonal_scaled,
        model.spectral_transform,
    )
    meridional_spectral = SpeedyWeather.SpeedyTransforms.transform(
        meridional_scaled,
        model.spectral_transform,
    )
    vorticity = SpeedyWeather.get_step(vars.prognostic.vorticity, 1)
    divergence = SpeedyWeather.get_step(vars.prognostic.divergence, 1)
    radius = model.geometry.radius[]
    SpeedyWeather.SpeedyTransforms.curl!(
        vorticity,
        zonal_spectral,
        meridional_spectral,
        model.spectral_transform;
        radius,
    )
    SpeedyWeather.SpeedyTransforms.divergence!(
        divergence,
        zonal_spectral,
        meridional_spectral,
        model.spectral_transform;
        radius,
    )
    SpeedyWeather.set!(
        vars,
        model;
        humidity = humidity_function,
        lf = 1,
        static_func = false,
    )
    return nothing
end

function _atmosphere_initial_conditions(config::ExperimentConfig, spectral_grid)
    if config.atmosphere_initial_conditions in (
        :era5_monthly,
        :era5_instantaneous,
    )
        return (; era5 = ERA5Atmosphere(config))
    end
    return SpeedyWeather.InitialConditions(spectral_grid, SpeedyWeather.PrimitiveWet)
end

function _atmosphere_initial_condition_provenance(initial_conditions, config)
    if config.atmosphere_initial_conditions == :analytic
        return (
            era5_source_mode = "not_applicable",
            era5_source_time = "not_applicable",
            era5_source_time_index = 0,
            era5_pressure_level_count = 0,
            era5_pressure_levels_hpa = "not_applicable",
            era5_pressure_levels_path = "not_applicable",
            era5_single_levels_path = "not_applicable",
            era5_pressure_levels_sha256 = "not_applicable",
            era5_single_levels_sha256 = "not_applicable",
            era5_vertical_interpolation = "not_applicable",
            era5_humidity_column_adjustment = "not_applicable",
        )
    end

    hasproperty(initial_conditions, :era5) || error(
        "configured ERA5 initial conditions are absent from the atmosphere model",
    )
    ic = initial_conditions.era5
    return (
        era5_source_mode = String(ic.source_mode),
        era5_source_time = string(ic.source_time),
        era5_source_time_index = ic.source_time_index,
        era5_pressure_level_count = length(ic.pressure_levels_hpa),
        era5_pressure_levels_hpa = join(string.(ic.pressure_levels_hpa), ","),
        era5_pressure_levels_path = ic.pressure_levels_path,
        era5_single_levels_path = ic.single_levels_path,
        era5_pressure_levels_sha256 = ic.pressure_levels_sha256,
        era5_single_levels_sha256 = ic.single_levels_sha256,
        era5_vertical_interpolation =
            "linear_in_log_pressure_temperature_and_wind_loglinear_specific_humidity",
        era5_humidity_column_adjustment =
            "column_scaled_to_ERA5_TCWV_with_98_percent_saturation_cap",
    )
end
