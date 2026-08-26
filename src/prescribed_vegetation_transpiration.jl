"""
Prescribed ERA5 vegetation evapotranspiration for the ReadyESM Terrarium stack.

The scheme deliberately does not activate Terrarium's unfinished prognostic
vegetation-carbon/canopy-water model. ERA5 high/low vegetation cover and LAI
are fixed boundary data. Bare/under-canopy evaporation remains a top-soil
flux, while transpiration is removed conservatively from the root zone in
proportion to static root fraction times current plant-available water.
"""
struct ERA5PrescribedVegetationEvapotranspiration{
    NF,
    GR <: Terrarium.AbstractGroundEvaporationResistanceFactor,
    PAW <: Terrarium.AbstractPlantAvailableWater,
} <: Terrarium.AbstractEvapotranspiration{NF}
    ground_resistance::GR
    plant_available_water::PAW
    maximum_leaf_conductance::NF
    minimum_leaf_conductance::NF
    light_half_saturation::NF
    vpd_scale::NF
    canopy_extinction_coefficient::NF
end

Adapt.@adapt_structure ERA5PrescribedVegetationEvapotranspiration

function ERA5PrescribedVegetationEvapotranspiration(
    ::Type{NF};
    maximum_leaf_conductance = NF(2e-3),
    minimum_leaf_conductance = NF(5e-5),
    light_half_saturation = NF(100),
    vpd_scale = NF(2000),
    canopy_extinction_coefficient = NF(0.5),
    ground_resistance = Terrarium.SoilMoistureResistanceFactor(NF),
    plant_available_water = Terrarium.FieldCapacityLimitedPAW(NF),
) where {NF}
    return ERA5PrescribedVegetationEvapotranspiration(
        ground_resistance,
        plant_available_water,
        NF(maximum_leaf_conductance),
        NF(minimum_leaf_conductance),
        NF(light_half_saturation),
        NF(vpd_scale),
        NF(canopy_extinction_coefficient),
    )
end

function _terrarium_evapotranspiration(config::ExperimentConfig, ::Type{NF}) where {NF}
    config.terrarium_evapotranspiration == :bare_ground &&
        return Terrarium.BareGroundEvaporation(NF)
    config.terrarium_evapotranspiration == :era5_prescribed_vegetation || error(
        "unsupported Terrarium evapotranspiration $(config.terrarium_evapotranspiration)",
    )
    return ERA5PrescribedVegetationEvapotranspiration(
        NF;
        maximum_leaf_conductance = config.terrarium_max_leaf_conductance_ms,
        minimum_leaf_conductance = config.terrarium_min_leaf_conductance_ms,
        light_half_saturation = config.terrarium_light_half_saturation_wm2,
        vpd_scale = config.terrarium_vpd_scale_pa,
        canopy_extinction_coefficient =
            config.terrarium_canopy_extinction_coefficient,
    )
end

function Terrarium.variables(
    ::ERA5PrescribedVegetationEvapotranspiration{NF},
) where {NF}
    # Reuse Terrarium's canonical declarations verbatim. In particular,
    # `skin_temperature` is also declared by the surface-energy process and
    # Terrarium rejects otherwise-equivalent variables whose metadata differ.
    bare_ground_variables = Terrarium.variables(Terrarium.BareGroundEvaporation(NF))
    ground_evaporation_conductance = bare_ground_variables[1]
    evaporation_ground = bare_ground_variables[2]
    skin_temperature = bare_ground_variables[3]
    return (
    ground_evaporation_conductance,
    Terrarium.auxiliary(
        :transpiration_conductance,
        Terrarium.XY();
        desc = "Vegetation transpiration vapor conductance",
    ),
    evaporation_ground,
    Terrarium.auxiliary(
        :transpiration,
        Terrarium.XY();
        desc = "Vegetation transpiration as liquid-water-depth flux",
    ),
    Terrarium.auxiliary(
        :ground_water_flux,
        Terrarium.XY();
        desc = "Ground evaporation as liquid-water-depth flux",
    ),
    Terrarium.auxiliary(
        :transpiration_water_flux,
        Terrarium.XY();
        desc = "Transpiration as liquid-water-depth flux",
    ),
    Terrarium.auxiliary(
        :root_water_availability,
        Terrarium.XYZ();
        desc = "Layer plant-available-water fraction",
    ),
    Terrarium.auxiliary(
        :root_uptake_weight_sum,
        Terrarium.XY();
        desc = "Root-weighted plant-available-water sum",
    ),
    Terrarium.input(
        :vegetation_fraction,
        Terrarium.XY();
        default = zero(NF),
        desc = "Prescribed ERA5 vegetated area fraction",
    ),
    Terrarium.input(
        :leaf_area_index,
        Terrarium.XY();
        default = zero(NF),
        desc = "Prescribed ERA5 LAI per vegetated area",
    ),
    Terrarium.input(
        :root_fraction,
        Terrarium.XYZ();
        default = zero(NF),
        desc = "Static fraction of roots in each soil layer",
    ),
    skin_temperature,
    )
end

Base.@propagate_inbounds function _era5_prescribed_surface_humidity_fluxes(
    i,
    j,
    grid,
    fields,
    evapotranspiration::ERA5PrescribedVegetationEvapotranspiration{NF},
    constants::Terrarium.PhysicalConstants,
    atmosphere::Terrarium.AbstractAtmosphere,
) where {NF}
    skin_temperature = fields.skin_temperature[i, j]
    humidity_difference = Terrarium.compute_specific_humidity_difference(
        i,
        j,
        grid,
        fields,
        atmosphere,
        constants,
        skin_temperature,
    )
    ground_humidity_flux = Terrarium.humidity_flux(
        evapotranspiration,
        humidity_difference,
        fields.ground_evaporation_conductance[i, j],
    )
    transpiration_humidity_flux = max(
        Terrarium.humidity_flux(
            evapotranspiration,
            humidity_difference,
            fields.transpiration_conductance[i, j],
        ),
        zero(NF),
    )
    return ground_humidity_flux, transpiration_humidity_flux
end

Base.@propagate_inbounds function Terrarium.compute_surface_humidity_flux(
    i,
    j,
    grid,
    fields,
    evapotranspiration::ERA5PrescribedVegetationEvapotranspiration,
    constants::Terrarium.PhysicalConstants,
    atmosphere::Terrarium.AbstractAtmosphere,
    args...,
)
    ground_flux, transpiration_flux = _era5_prescribed_surface_humidity_fluxes(
        i,
        j,
        grid,
        fields,
        evapotranspiration,
        constants,
        atmosphere,
    )
    return ground_flux + transpiration_flux
end

Base.@propagate_inbounds function Terrarium.ground_evapotranspiration_flux(
    i,
    j,
    grid,
    fields,
    ::ERA5PrescribedVegetationEvapotranspiration,
    args...,
)
    return fields.ground_water_flux[i, j] + fields.transpiration_water_flux[i, j]
end

KernelAbstractions.@kernel function _era5_root_water_availability_kernel!(
    root_water_availability,
    grid,
    fields,
    plant_available_water,
    stratigraphy,
    hydrology,
    biogeochemistry,
)
    i, j, k = @index(Global, NTuple)
    value = Terrarium.compute_plant_available_water(
        i,
        j,
        k,
        grid,
        fields,
        plant_available_water,
        stratigraphy,
        hydrology,
        biogeochemistry,
    )
    @inbounds root_water_availability[i, j, k] = clamp(value, zero(value), one(value))
end

KernelAbstractions.@kernel function _era5_prescribed_vegetation_flux_kernel!(
    outputs,
    grid,
    fields,
    evapotranspiration,
    constants,
    atmosphere,
    soil,
    snow,
)
    i, j = @index(Global, NTuple)
    NF = eltype(grid)
    @inbounds begin
        vegetation_fraction = clamp(fields.vegetation_fraction[i, j], zero(NF), one(NF))
        leaf_area_index = clamp(fields.leaf_area_index[i, j], zero(NF), NF(10))
        root_uptake_weight_sum = zero(NF)
        for k in 1:Terrarium.get_field_grid(grid).Nz
            root_uptake_weight_sum +=
                max(fields.root_fraction[i, j, k], zero(NF)) *
                clamp(fields.root_water_availability[i, j, k], zero(NF), one(NF))
        end

        aerodynamic_resistance = Terrarium.aerodynamic_resistance(
            i,
            j,
            grid,
            fields,
            atmosphere,
        )
        ground_resistance = Terrarium.ground_evaporation_resistance_factor(
            i,
            j,
            grid,
            fields,
            evapotranspiration.ground_resistance,
            soil,
        )
        shortwave = max(
            Terrarium.shortwave_down(i, j, grid, fields, atmosphere),
            zero(NF),
        )
        light_factor = shortwave /
            (shortwave + evapotranspiration.light_half_saturation)
        vpd = max(
            Terrarium.compute_vapor_pressure_deficit(
                i,
                j,
                grid,
                fields,
                atmosphere,
                constants,
            ),
            zero(NF),
        )
        vpd_factor = evapotranspiration.vpd_scale /
            (evapotranspiration.vpd_scale + vpd)
        leaf_conductance = evapotranspiration.minimum_leaf_conductance +
            (evapotranspiration.maximum_leaf_conductance -
             evapotranspiration.minimum_leaf_conductance) *
            light_factor * vpd_factor
        canopy_conductance = leaf_conductance * leaf_area_index *
            clamp(root_uptake_weight_sum, zero(NF), one(NF))
        transpiration_conductance = ifelse(
            canopy_conductance > eps(NF),
            vegetation_fraction /
            (aerodynamic_resistance + inv(max(canopy_conductance, eps(NF)))),
            zero(NF),
        )
        ground_exposure = one(NF) - vegetation_fraction *
            (one(NF) - exp(
                -evapotranspiration.canopy_extinction_coefficient *
                leaf_area_index,
            ))
        ground_evaporation_conductance =
            ground_exposure * ground_resistance / aerodynamic_resistance

        outputs.ground_evaporation_conductance[i, j, 1] =
            ground_evaporation_conductance
        outputs.transpiration_conductance[i, j, 1] = transpiration_conductance
        outputs.root_uptake_weight_sum[i, j, 1] = root_uptake_weight_sum

        conductance_fields = merge(
            fields,
            (;
                ground_evaporation_conductance =
                    outputs.ground_evaporation_conductance,
                transpiration_conductance = outputs.transpiration_conductance,
            ),
        )
        Terrarium.compute_evapotranspiration_fluxes!(
            outputs,
            i,
            j,
            grid,
            conductance_fields,
            evapotranspiration,
            constants,
            atmosphere,
            snow,
        )
    end
end

Base.@propagate_inbounds function Terrarium.compute_evapotranspiration_fluxes!(
    outputs,
    i,
    j,
    grid,
    fields,
    evapotranspiration::ERA5PrescribedVegetationEvapotranspiration{NF},
    constants::Terrarium.PhysicalConstants,
    atmosphere::Terrarium.AbstractAtmosphere,
    snow::Terrarium.Optional{Terrarium.AbstractSnow} = nothing,
    args...,
) where {NF}
    ground_humidity_flux, transpiration_humidity_flux =
        _era5_prescribed_surface_humidity_fluxes(
            i,
            j,
            grid,
            fields,
            evapotranspiration,
            constants,
            atmosphere,
        )
    snow_free_fraction = one(NF) -
        Terrarium.snow_cover_fraction(i, j, grid, fields, snow)
    density_ratio =
        Terrarium.air_density(i, j, grid, fields, atmosphere, constants) /
        constants.material.density_water
    ground_water_flux =
        snow_free_fraction * density_ratio * ground_humidity_flux
    transpiration_water_flux =
        snow_free_fraction * density_ratio * transpiration_humidity_flux

    outputs.evaporation_ground[i, j, 1] = ground_water_flux
    outputs.transpiration[i, j, 1] = transpiration_water_flux
    outputs.ground_water_flux[i, j, 1] = ground_water_flux
    outputs.transpiration_water_flux[i, j, 1] = transpiration_water_flux
    return outputs
end

function Terrarium.compute_auxiliary!(
    state,
    grid,
    evapotranspiration::ERA5PrescribedVegetationEvapotranspiration,
    ::Terrarium.NoCanopyInterception,
    constants::Terrarium.PhysicalConstants,
    atmosphere::Terrarium.AbstractAtmosphere,
    soil::Terrarium.AbstractSoil,
    snow::Terrarium.Optional{Terrarium.AbstractSnow},
    args...,
)
    hydrology = Terrarium.get_hydrology(soil)
    stratigraphy = Terrarium.get_stratigraphy(soil)
    biogeochemistry = Terrarium.get_biogeochemistry(soil)
    availability_output = state.root_water_availability
    availability_fields = Terrarium.get_fields(
        state,
        evapotranspiration,
        soil;
        except = (; root_water_availability = availability_output),
    )
    Terrarium.launch!(
        grid,
        Terrarium.XYZ,
        _era5_root_water_availability_kernel!,
        availability_output,
        availability_fields,
        evapotranspiration.plant_available_water,
        stratigraphy,
        hydrology,
        biogeochemistry,
    )

    outputs = (
        ground_evaporation_conductance = state.ground_evaporation_conductance,
        transpiration_conductance = state.transpiration_conductance,
        evaporation_ground = state.evaporation_ground,
        transpiration = state.transpiration,
        ground_water_flux = state.ground_water_flux,
        transpiration_water_flux = state.transpiration_water_flux,
        root_uptake_weight_sum = state.root_uptake_weight_sum,
    )
    fields = Terrarium.get_fields(
        state,
        evapotranspiration,
        atmosphere,
        soil,
        snow;
        except = outputs,
    )
    Terrarium.launch!(
        grid,
        Terrarium.XY,
        _era5_prescribed_vegetation_flux_kernel!,
        outputs,
        fields,
        evapotranspiration,
        constants,
        atmosphere,
        soil,
        snow,
    )
    return nothing
end

Base.@propagate_inbounds function Terrarium.forcing(
    i,
    j,
    k,
    grid,
    clock,
    fields,
    ::ERA5PrescribedVegetationEvapotranspiration,
    ::Terrarium.SoilHydrology{NF, Terrarium.RichardsEq},
    ::Terrarium.PhysicalConstants,
) where {NF}
    field_grid = Terrarium.get_field_grid(grid)
    layer_thickness = Terrarium.Δzᵃᵃᶜ(i, j, k, field_grid)
    ground_flux = fields.ground_water_flux[i, j]
    transpiration_flux = fields.transpiration_water_flux[i, j]
    uptake_weight_sum = fields.root_uptake_weight_sum[i, j]
    uptake_fraction = ifelse(
        uptake_weight_sum > eps(NF),
        max(fields.root_fraction[i, j, k], zero(NF)) *
        clamp(fields.root_water_availability[i, j, k], zero(NF), one(NF)) /
        max(uptake_weight_sum, eps(NF)),
        zero(NF),
    )
    ground_tendency = -ground_flux / layer_thickness * (k == field_grid.Nz)
    transpiration_tendency = -transpiration_flux * uptake_fraction / layer_thickness
    return ground_tendency + transpiration_tendency
end
