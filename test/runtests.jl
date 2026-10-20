using ReadyESM
using Test

include("ordered_ocean_grid.jl")
include("standalone_atmosphere_restart.jl")
include("carbon_stage0a.jl")
include("carbon_stage0b.jl")
include("carbon_stage0c.jl")
include("carbon_emissions.jl")
include("carbon_diagnostics.jl")
include("carbon_inventory_geometry.jl")

config = load_config(joinpath(ReadyESM.PROJECT_ROOT, "config", "production.yml");
    check_input_files = false)

@test config.device == :gpu
@test (config.truncation, config.nlayers) == (31, 27)
@test config.ocean_grid == :tripolar_1degree
@test (config.ocean_nlongitude, config.ocean_nlatitude, config.ocean_nlayers) ==
      (360, 180, 60)
@test config.ocean_dynamics
@test config.sea_ice_dynamics
@test config.land_model == :terrarium
@test config.terrarium_precision == :float64
@test config.terrarium_soil_layers == 16
@test config.atmosphere_vertical_diffusion ==
      :speedyweather_0_21_1_zero_operator_control
@test config.atmosphere_large_scale_precipitation == :net_column_conservative
@test config.atmosphere_saturation_law == :mixed_phase
@test config.atmosphere_convection == :betts_miller_constant_rh
@test config.atmosphere_convection_timescale_seconds == 7200
@test config.atmosphere_initial_spectral_truncation
@test !ExperimentConfig().atmosphere_initial_spectral_truncation
include("convection_timescale.jl")
include("finite_multiplume_convection.jl")
include("convective_momentum_transport.jl")
include("optional_land.jl")
@test config.forcing.radiation == :rrtmgp_all_sky
@test config.forcing.co2_ppm == 420
@test config.forcing.cloud_liquid_water_path_gm2 == 75
@test config.forcing.cloud_ice_water_path_gm2 == 25
@test config.forcing.aerosol_optical_depth_550nm == 0

diagnostic_source = read(joinpath(ReadyESM.PROJECT_ROOT, "src", "diagnostics.jl"), String)
dynamic_source = read(joinpath(ReadyESM.PROJECT_ROOT, "src", "dynamic_ocean.jl"), String)
for contract in (
    "RRTMGP.direct_sw_surface_albedo",
    "RRTMGP.diffuse_sw_surface_albedo",
    "RRTMGP direct and diffuse surface albedos differ",
    "RRTMGP surface albedo is not band invariant",
)
    @test occursin(contract, diagnostic_source)
end
for contract in (
    "atmosphere_radiative_surface_albedo",
    "retained_rrtmgp_direct_diffuse_band_invariant_solver_state",
    "radiative surface albedo is outside [0, 1]",
)
    @test occursin(contract, dynamic_source)
end

et = ReadyESM.ERA5PrescribedVegetationEvapotranspiration(Float32)
land_dispatch = which(
    ReadyESM.Terrarium.compute_auxiliary!,
    Tuple{
        Any,
        Any,
        typeof(et),
        typeof(ReadyESM.Terrarium.NoCanopyInterception(Float32)),
        ReadyESM.Terrarium.PhysicalConstants{Float32},
        ReadyESM.Terrarium.AbstractAtmosphere,
        ReadyESM.Terrarium.AbstractSoil,
        Nothing,
        Nothing,
    },
)
@test land_dispatch.module === ReadyESM

include("config_input_files.jl")
include("calendar_centred_ecco.jl")
include("darcy_boundary.jl")
include("sampled_global_diagnostics.jl")
include("sampled_history_output.jl")
include("ocean_history_bounds.jl")
include("native_runoff.jl")
include("bs08_runoff_enthalpy.jl")
include("tripolar_runoff_delivery.jl")
include("precipitation_diagnostics.jl")
include("precipitation_pressure.jl")
include("sea_ice_branch_types.jl")
include("sea_ice_thermal_conservation.jl")
include("land_precision.jl")
include("mixed_phase_saturation.jl")
include("atmosphere_spectral_truncation.jl")
