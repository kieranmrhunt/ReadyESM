using ReadyESM
using Test

config = load_config(joinpath(ReadyESM.PROJECT_ROOT, "config", "production.yml"))

@test config.device == :gpu
@test (config.truncation, config.nlayers) == (31, 27)
@test config.ocean_grid == :tripolar_1degree
@test (config.ocean_nlongitude, config.ocean_nlatitude, config.ocean_nlayers) ==
      (360, 180, 60)
@test config.ocean_dynamics
@test config.sea_ice_dynamics
@test config.land_model == :terrarium
@test config.terrarium_soil_layers == 16
@test config.atmosphere_vertical_diffusion ==
      :speedyweather_0_21_1_zero_operator_control
@test config.atmosphere_large_scale_precipitation == :upstream_observed
@test config.atmosphere_convection == :betts_miller_constant_rh
@test config.forcing.radiation == :rrtmgp_all_sky
@test config.forcing.co2_ppm == 420
@test config.forcing.aerosol_optical_depth_550nm == 0
