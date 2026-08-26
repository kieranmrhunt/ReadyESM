const ERA5_LAND_LAYER_INTERFACES_M = Float32[0, 0.07, 0.28, 1.0, 2.89]

# Saturated volumetric water contents used by the ERA5 IFS/H-TESSEL soil
# classes (1: coarse, 2: medium, 3: medium-fine, 4: fine, 5: very fine,
# 6: extended tropical organic, 7: tropical organic).  Converting ERA5 VWC
# to relative saturation with its *source* soil class preserves the reanalysis
# wetness state without pretending Terrarium's homogeneous soil has the same
# water-holding capacity everywhere.
const ERA5_SOIL_SATURATED_VWC = Float32[
    0.403,
    0.439,
    0.430,
    0.520,
    0.614,
    0.766,
    0.439,
]

struct ERA5LandState{V, A3, A2}
    latitude::V
    longitude::V
    volumetric_soil_water::A3
    soil_temperature::A3
    land_fraction::A2
    soil_type::A2
    skin_temperature::A2
    high_vegetation_cover::A2
    low_vegetation_cover::A2
    high_vegetation_lai::A2
    low_vegetation_lai::A2
    source_time::DateTime
    source_time_index::Int
    path::String
    sha256::String
end

function _era5_exact_selected_time(dataset, expected_time, description)
    _, times = _era5_time_coordinate(dataset)
    matches = findall(==(expected_time), times)
    length(matches) == 1 || error(
        "$description must contain model start timestamp $expected_time " *
        "exactly once; found $(length(matches)) matches",
    )
    index = only(matches)
    return index, times[index]
end

function _era5_land_matrix(dataset, name, time_index)
    return permutedims(
        _era5_array(dataset[name], :, :, time_index),
        (2, 1),
    )
end

function ERA5LandState(config::ExperimentConfig)
    config.terrarium_initial_conditions == :era5_instantaneous || error(
        "ERA5LandState requires terrarium_initial_conditions=era5_instantaneous",
    )

    latitude = Float32[]
    longitude = Float32[]
    soil_water = Array{Float32}(undef, 0, 0, 0)
    soil_temperature = similar(soil_water)
    land_fraction = Array{Float32}(undef, 0, 0)
    soil_type = similar(land_fraction)
    skin_temperature = similar(land_fraction)
    high_vegetation_cover = similar(land_fraction)
    low_vegetation_cover = similar(land_fraction)
    high_vegetation_lai = similar(land_fraction)
    low_vegetation_lai = similar(land_fraction)
    source_time = DateTime(0)
    source_time_index = 0

    NCDataset(config.era5_land_state_path, "r") do dataset
        source_time_index, source_time = _era5_exact_selected_time(
            dataset,
            config.start_date,
            "ERA5 instantaneous land-state file",
        )
        _era5_require_units(dataset["latitude"], "land latitude", ("degreesnorth",))
        _era5_require_units(dataset["longitude"], "land longitude", ("degreeseast",))
        _era5_require_units(dataset["lsm"], "land-sea mask", ("(0-1)", "1"))
        _era5_require_units(dataset["slt"], "soil type", ("(codetable4.213)",))
        _era5_require_units(dataset["skt"], "skin temperature", ("k",))
        for name in ("swvl1", "swvl2", "swvl3", "swvl4")
            _era5_require_units(
                dataset[name],
                name,
                ("m^3m^-3", "m3m-3"),
            )
        end
        for name in ("stl1", "stl2", "stl3", "stl4")
            _era5_require_units(dataset[name], name, ("k",))
        end
        for name in ("cvh", "cvl")
            _era5_require_units(dataset[name], name, ("(0-1)", "1"))
        end
        for name in ("lai_hv", "lai_lv")
            _era5_require_units(dataset[name], name, ("m^2m^-2", "m2m-2"))
        end
        for name in (
            "lsm", "slt", "skt",
            "swvl1", "swvl2", "swvl3", "swvl4",
            "stl1", "stl2", "stl3", "stl4",
            "cvh", "cvl", "lai_hv", "lai_lv",
        )
            _era5_require_instantaneous(dataset[name], "land-state $name")
        end

        latitude = _era5_array(dataset["latitude"], :)
        longitude = _era5_array(dataset["longitude"], :)
        water_layers = ntuple(
            layer -> _era5_land_matrix(dataset, "swvl$layer", source_time_index),
            4,
        )
        temperature_layers = ntuple(
            layer -> _era5_land_matrix(dataset, "stl$layer", source_time_index),
            4,
        )
        soil_water = cat(water_layers...; dims = 3)
        soil_water = permutedims(soil_water, (3, 1, 2))
        soil_temperature = cat(temperature_layers...; dims = 3)
        soil_temperature = permutedims(soil_temperature, (3, 1, 2))
        land_fraction = _era5_land_matrix(dataset, "lsm", source_time_index)
        soil_type = _era5_land_matrix(dataset, "slt", source_time_index)
        skin_temperature = _era5_land_matrix(dataset, "skt", source_time_index)
        high_vegetation_cover = _era5_land_matrix(dataset, "cvh", source_time_index)
        low_vegetation_cover = _era5_land_matrix(dataset, "cvl", source_time_index)
        high_vegetation_lai = _era5_land_matrix(dataset, "lai_hv", source_time_index)
        low_vegetation_lai = _era5_land_matrix(dataset, "lai_lv", source_time_index)
    end

    all(diff(longitude) .> 0) || error("ERA5 land longitude must increase")
    all(diff(latitude) .< 0) || error("ERA5 land latitude must decrease")
    expected_2d = (length(latitude), length(longitude))
    expected_3d = (4, expected_2d...)
    size(soil_water) == expected_3d || error(
        "unexpected ERA5 soil-water shape $(size(soil_water)); expected $expected_3d",
    )
    size(soil_temperature) == expected_3d || error(
        "unexpected ERA5 soil-temperature shape $(size(soil_temperature)); expected $expected_3d",
    )
    all(size(field) == expected_2d for field in (
        land_fraction,
        soil_type,
        skin_temperature,
        high_vegetation_cover,
        low_vegetation_cover,
        high_vegetation_lai,
        low_vegetation_lai,
    )) || error("ERA5 land-state variables do not share one coordinate shape")

    _era5_validate_finite_range(land_fraction, "land-sea mask", 0, 1)
    _era5_validate_finite_range(soil_type, "soil type", 0, 7)
    all(value -> isinteger(value), soil_type) ||
        error("ERA5 soil type contains non-integer category codes")
    _era5_validate_finite_range(soil_water, "volumetric soil water", -0.005, 0.8)
    _era5_validate_finite_range(soil_temperature, "soil temperature", 180, 340)
    _era5_validate_finite_range(skin_temperature, "land skin temperature", 180, 340)
    _era5_validate_finite_range(high_vegetation_cover, "high vegetation cover", 0, 1)
    _era5_validate_finite_range(low_vegetation_cover, "low vegetation cover", 0, 1)
    _era5_validate_finite_range(high_vegetation_lai, "high vegetation LAI", 0, 10)
    _era5_validate_finite_range(low_vegetation_lai, "low vegetation LAI", 0, 10)

    return ERA5LandState(
        latitude,
        longitude,
        soil_water,
        soil_temperature,
        land_fraction,
        soil_type,
        skin_temperature,
        high_vegetation_cover,
        low_vegetation_cover,
        high_vegetation_lai,
        low_vegetation_lai,
        source_time,
        source_time_index,
        config.era5_land_state_path,
        _era5_sha256(config.era5_land_state_path),
    )
end

function _era5_land_stencil(source::ERA5LandState, longitude, latitude)
    i0, i1, fi, j0, j1, fj =
        _era5_horizontal_indices(source, longitude, latitude)
    return (
        ((j0, i0), (j0, i1), (j1, i0), (j1, i1)),
        ((1 - fi) * (1 - fj), fi * (1 - fj), (1 - fi) * fj, fi * fj),
    )
end

function _era5_nearest_land_value(field, source::ERA5LandState, longitude, latitude)
    i_center = argmin(abs.(source.longitude .- mod(longitude, 360)))
    j_center = argmin(abs.(source.latitude .- latitude))
    nlongitude = length(source.longitude)
    nlatitude = length(source.latitude)
    for radius in 0:12
        best_distance = Inf
        best_value = NaN
        for dj in -radius:radius, di in -radius:radius
            max(abs(di), abs(dj)) == radius || continue
            i = mod1(i_center + di, nlongitude)
            j = clamp(j_center + dj, 1, nlatitude)
            source.land_fraction[j, i] > 1e-6 || continue
            value = field[j, i]
            isfinite(value) || continue
            distance = di^2 + dj^2
            if distance < best_distance
                best_distance = distance
                best_value = value
            end
        end
        isfinite(best_value) && return best_value
    end
    error("no valid ERA5 land value within 12 degrees of ($longitude, $latitude)")
end

function _era5_land_bilinear(field, source::ERA5LandState, longitude, latitude)
    indices, weights = _era5_land_stencil(source, longitude, latitude)
    numerator = 0.0
    denominator = 0.0
    for ((j, i), weight) in zip(indices, weights)
        land_weight = weight * clamp(source.land_fraction[j, i], 0, 1)
        value = field[j, i]
        if isfinite(value) && land_weight > 0
            numerator += land_weight * value
            denominator += land_weight
        end
    end
    denominator > 0 && return numerator / denominator
    return _era5_nearest_land_value(field, source, longitude, latitude)
end

function _era5_land_nearest_category(field, source::ERA5LandState, longitude, latitude)
    indices, weights = _era5_land_stencil(source, longitude, latitude)
    scores = map(indices, weights) do (j, i), weight
        value = field[j, i]
        valid = isfinite(value) && source.land_fraction[j, i] > 1e-6
        return valid ? weight * source.land_fraction[j, i] : -Inf
    end
    maximum(scores) > -Inf && return field[indices[argmax(scores)]...]
    return _era5_nearest_land_value(field, source, longitude, latitude)
end

function _piecewise_constant_integral(values, interfaces, shallow_depth, deep_depth)
    integral = 0.0
    for layer in eachindex(values)
        left = max(shallow_depth, interfaces[layer])
        right = min(deep_depth, interfaces[layer + 1])
        right > left && (integral += (right - left) * values[layer])
    end
    if deep_depth > interfaces[end]
        left = max(shallow_depth, interfaces[end])
        deep_depth > left && (integral += (deep_depth - left) * values[end])
    end
    return integral
end

function _era5_relative_saturation_integral(
    layer_saturation,
    shallow_depth,
    deep_depth,
    water_table_depth,
)
    interfaces = ERA5_LAND_LAYER_INTERFACES_M
    source_bottom = Float64(interfaces[end])
    upper_bottom = min(deep_depth, source_bottom)
    integral = shallow_depth < upper_bottom ?
        _piecewise_constant_integral(
            layer_saturation,
            interfaces,
            shallow_depth,
            upper_bottom,
        ) : 0.0

    transition_top = max(shallow_depth, source_bottom)
    transition_bottom = min(deep_depth, water_table_depth)
    if transition_bottom > transition_top
        deepest_saturation = layer_saturation[end]
        slope = (1 - deepest_saturation) / (water_table_depth - source_bottom)
        primitive(depth) =
            deepest_saturation * depth + slope * (depth - source_bottom)^2 / 2
        integral += primitive(transition_bottom) - primitive(transition_top)
    end
    saturated_top = max(shallow_depth, water_table_depth)
    deep_depth > saturated_top && (integral += deep_depth - saturated_top)
    return integral
end

function _terrarium_vertical_face_depths(column_grid)
    field_grid = Terrarium.get_field_grid(column_grid)
    faces = Oceananigans.Grids.znodes(
        field_grid,
        Oceananigans.Center(),
        Oceananigans.Center(),
        Oceananigans.Face(),
    )
    host_faces = Oceananigans.Architectures.on_architecture(
        Oceananigans.CPU(),
        faces,
    )
    return -Float64.(collect(host_faces))
end

function _static_exponential_root_fractions(face_depths, a, b)
    a > 0 || throw(ArgumentError("root-distribution parameter a must be positive"))
    b > 0 || throw(ArgumentError("root-distribution parameter b must be positive"))
    fractions = Float64[]
    for layer in 1:(length(face_depths) - 1)
        shallow = min(face_depths[layer], face_depths[layer + 1])
        deep = max(face_depths[layer], face_depths[layer + 1])
        push!(fractions, 0.5 * (
            exp(-a * shallow) - exp(-a * deep) +
            exp(-b * shallow) - exp(-b * deep)
        ))
    end
    total = sum(fractions)
    total > 0 || error("root distribution has zero integral over the soil column")
    fractions ./= total
    all(isfinite, fractions) || error("root distribution is non-finite")
    all(>=(0), fractions) || error("root distribution contains a negative fraction")
    return fractions
end

function _era5_land_initializers(
    config::ExperimentConfig,
    spectral_grid,
    column_grid,
    water_table_elevation,
)
    source = ERA5LandState(config)
    longitudes, latitudes =
        SpeedyWeather.RingGrids.get_londlatds(spectral_grid.grid)
    mask = Bool.(Array(column_grid.mask.data))
    length(mask) == length(longitudes) == length(latitudes) || error(
        "Terrarium mask and SpeedyWeather coordinates do not share one grid",
    )
    land_points = findall(mask)
    nland = length(land_points)
    face_depths = _terrarium_vertical_face_depths(column_grid)
    nlayer = length(face_depths) - 1
    water_table_depth = -Float64(water_table_elevation)
    water_table_depth > ERA5_LAND_LAYER_INTERFACES_M[end] || error(
        "Terrarium water table must lie below ERA5's 2.89 m soil column",
    )

    temperature = Array{Float32}(undef, nland, 1, nlayer)
    saturation = similar(temperature)
    skin_temperature = Array{Float32}(undef, nland, 1, 1)
    prescribed_vegetation =
        config.terrarium_evapotranspiration == :era5_prescribed_vegetation
    vegetation_fraction = prescribed_vegetation ?
        Array{Float32}(undef, nland, 1, 1) : Array{Float32}(undef, 0, 1, 1)
    leaf_area_index = similar(vegetation_fraction)
    root_fraction = prescribed_vegetation ?
        Array{Float32}(undef, nland, 1, nlayer) :
        Array{Float32}(undef, 0, 1, nlayer)
    root_layer_fractions = prescribed_vegetation ?
        _static_exponential_root_fractions(
            face_depths,
            config.terrarium_root_a_m1,
            config.terrarium_root_b_m1,
        ) : Float64[]
    zero_soil_type_count = 0
    source_saturation_clip_count = 0

    for (column, grid_point) in enumerate(land_points)
        longitude = longitudes[grid_point]
        latitude = latitudes[grid_point]
        soil_type_raw = _era5_land_nearest_category(
            source.soil_type,
            source,
            longitude,
            latitude,
        )
        soil_type = round(Int, soil_type_raw)
        if soil_type in eachindex(ERA5_SOIL_SATURATED_VWC)
            saturated_vwc = Float64(ERA5_SOIL_SATURATED_VWC[soil_type])
        else
            zero_soil_type_count += 1
            saturated_vwc = 0.49
        end
        source_water = Float64[
            _era5_land_bilinear(
                view(source.volumetric_soil_water, layer, :, :),
                source,
                longitude,
                latitude,
            ) for layer in 1:4
        ]
        raw_saturation = source_water ./ saturated_vwc
        source_saturation_clip_count += count(value -> value < 0.01 || value > 1, raw_saturation)
        layer_saturation = clamp.(raw_saturation, 0.01, 1.0)
        source_temperature = Float64[
            _era5_land_bilinear(
                view(source.soil_temperature, layer, :, :),
                source,
                longitude,
                latitude,
            ) for layer in 1:4
        ]

        for layer in 1:nlayer
            shallow_depth = min(face_depths[layer], face_depths[layer + 1])
            deep_depth = max(face_depths[layer], face_depths[layer + 1])
            thickness = deep_depth - shallow_depth
            saturation[column, 1, layer] = Float32(clamp(
                _era5_relative_saturation_integral(
                    layer_saturation,
                    shallow_depth,
                    deep_depth,
                    water_table_depth,
                ) / thickness,
                0.01,
                1.0,
            ))
            temperature[column, 1, layer] = Float32(
                _piecewise_constant_integral(
                    source_temperature,
                    ERA5_LAND_LAYER_INTERFACES_M,
                    shallow_depth,
                    deep_depth,
                ) / thickness - 273.15,
            )
        end
        skin_temperature[column, 1, 1] = Float32(
            _era5_land_bilinear(
                source.skin_temperature,
                source,
                longitude,
                latitude,
            ) - 273.15,
        )
        if prescribed_vegetation
            high_cover = clamp(
                _era5_land_bilinear(
                    source.high_vegetation_cover,
                    source,
                    longitude,
                    latitude,
                ),
                0,
                1,
            )
            low_cover = clamp(
                _era5_land_bilinear(
                    source.low_vegetation_cover,
                    source,
                    longitude,
                    latitude,
                ),
                0,
                1,
            )
            total_cover = clamp(high_cover + low_cover, 0, 1)
            grid_area_lai = high_cover * _era5_land_bilinear(
                source.high_vegetation_lai,
                source,
                longitude,
                latitude,
            ) + low_cover * _era5_land_bilinear(
                source.low_vegetation_lai,
                source,
                longitude,
                latitude,
            )
            vegetation_fraction[column, 1, 1] = Float32(total_cover)
            leaf_area_index[column, 1, 1] = Float32(clamp(
                total_cover > eps(Float64) ? grid_area_lai / total_cover : 0,
                0,
                10,
            ))
            for layer in 1:nlayer
                root_fraction[column, 1, layer] =
                    Float32(root_layer_fractions[layer])
            end
        end
    end

    all(isfinite, temperature) || error("remapped ERA5 soil temperature is non-finite")
    all(isfinite, saturation) || error("remapped ERA5 soil saturation is non-finite")
    all(isfinite, skin_temperature) || error("remapped ERA5 skin temperature is non-finite")
    all((-100 .<= temperature) .& (temperature .<= 80)) ||
        error("remapped ERA5 soil temperature is implausible: $(extrema(temperature))")
    all((0.01f0 .<= saturation) .& (saturation .<= 1f0)) ||
        error("remapped ERA5 soil saturation is out of bounds: $(extrema(saturation))")
    all((-100 .<= skin_temperature) .& (skin_temperature .<= 80)) ||
        error("remapped ERA5 skin temperature is implausible: $(extrema(skin_temperature))")
    if prescribed_vegetation
        all((0f0 .<= vegetation_fraction) .& (vegetation_fraction .<= 1f0)) ||
            error("remapped ERA5 vegetation fraction is out of bounds")
        all((0f0 .<= leaf_area_index) .& (leaf_area_index .<= 10f0)) ||
            error("remapped ERA5 LAI is out of bounds")
        all(abs.(vec(sum(root_fraction; dims = 3)) .- 1) .<= 5f-6) ||
            error("remapped static root fractions do not sum to one")
    end

    architecture = Terrarium.architecture(column_grid)
    initializers = (;
        temperature = Terrarium.on_architecture(architecture, temperature),
        saturation_water_ice = Terrarium.on_architecture(architecture, saturation),
        skin_temperature = Terrarium.on_architecture(architecture, skin_temperature),
    )
    if prescribed_vegetation
        initializers = merge(initializers, (;
            vegetation_fraction = Terrarium.on_architecture(
                architecture,
                vegetation_fraction,
            ),
            leaf_area_index = Terrarium.on_architecture(
                architecture,
                leaf_area_index,
            ),
            root_fraction = Terrarium.on_architecture(
                architecture,
                root_fraction,
            ),
        ))
    end
    audit = (;
        land_column_count = nland,
        zero_soil_type_count,
        source_saturation_clip_count,
        temperature_extrema_celsius = extrema(temperature),
        saturation_extrema = extrema(saturation),
        skin_temperature_extrema_celsius = extrema(skin_temperature),
        vegetation_fraction_extrema = prescribed_vegetation ?
            extrema(vegetation_fraction) : (0f0, 0f0),
        leaf_area_index_extrema = prescribed_vegetation ?
            extrema(leaf_area_index) : (0f0, 0f0),
        root_fraction_sum_extrema = prescribed_vegetation ?
            extrema(vec(sum(root_fraction; dims = 3))) : (0f0, 0f0),
    )
    return (; initializers, source, audit)
end

function _land_initial_condition_provenance(config::ExperimentConfig)
    common = (;
        terrarium_initial_conditions = String(config.terrarium_initial_conditions),
        era5_land_state_path = "",
        era5_land_state_sha256 = "",
        era5_land_source_time = "",
        era5_land_source_time_index = 0,
        era5_land_layer_interfaces_m = join(ERA5_LAND_LAYER_INTERFACES_M, ","),
        era5_land_soil_moisture_conversion = "none",
        era5_land_deep_profile = "none",
        terrarium_evapotranspiration = String(config.terrarium_evapotranspiration),
        era5_vegetation_source = "none",
        era5_vegetation_cover_lai_combination = "none",
        terrarium_root_distribution = "none",
        terrarium_max_leaf_conductance_ms = NaN,
        terrarium_min_leaf_conductance_ms = NaN,
        terrarium_light_half_saturation_wm2 = NaN,
        terrarium_vpd_scale_pa = NaN,
        terrarium_canopy_extinction_coefficient = NaN,
    )
    config.land_model == :terrarium || return common
    config.terrarium_initial_conditions == :era5_instantaneous || return common
    source = ERA5LandState(config)
    land_state = merge(common, (;
        era5_land_state_path = source.path,
        era5_land_state_sha256 = source.sha256,
        era5_land_source_time = string(source.source_time),
        era5_land_source_time_index = source.source_time_index,
        era5_land_soil_moisture_conversion =
            "ERA5-H-TESSEL-soil-class-relative-saturation",
        era5_land_deep_profile =
            "layer4-linear-to-aligned-5m-water-table-then-saturated",
    ))
    config.terrarium_evapotranspiration == :era5_prescribed_vegetation ||
        return land_state
    return merge(land_state, (;
        era5_vegetation_source = "ERA5-cvh-cvl-lai_hv-lai_lv-instantaneous",
        era5_vegetation_cover_lai_combination =
            "fveg=cvh+cvl;LAI=(cvh*lai_hv+cvl*lai_lv)/fveg",
        terrarium_root_distribution =
            "normalized-0.5*(a*exp(-a*z)+b*exp(-b*z));" *
            "a=$(config.terrarium_root_a_m1)m-1;b=$(config.terrarium_root_b_m1)m-1",
        terrarium_max_leaf_conductance_ms =
            config.terrarium_max_leaf_conductance_ms,
        terrarium_min_leaf_conductance_ms =
            config.terrarium_min_leaf_conductance_ms,
        terrarium_light_half_saturation_wm2 =
            config.terrarium_light_half_saturation_wm2,
        terrarium_vpd_scale_pa = config.terrarium_vpd_scale_pa,
        terrarium_canopy_extinction_coefficient =
            config.terrarium_canopy_extinction_coefficient,
    ))
end
