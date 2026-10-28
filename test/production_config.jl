using ReadyESM, Test, Dates

@testset "Released production configuration" begin
    config = load_config(joinpath(ReadyESM.PROJECT_ROOT, "config", "production.yml");
        check_input_files = false)
    @test config.device == :gpu
    @test (config.truncation, config.nlayers) == (31, 27)
    @test config.ocean_grid == :tripolar_1degree
    @test (config.ocean_nlongitude, config.ocean_nlatitude, config.ocean_nlayers) == (360, 180, 60)
    @test config.ocean_dynamics
    @test config.sea_ice_dynamics
    @test config.land_model == :terrarium
    @test config.terrarium_precision == :float64
    @test config.terrarium_soil_layers == 16
    @test config.atmosphere_vertical_diffusion == :holtslag_boville_k_profile_v1
    @test config.atmosphere_pbl_diffusivity_scale == 1.0
    @test config.atmosphere_pbl_mix_momentum
    @test !config.atmosphere_pbl_mix_heat
    @test !config.atmosphere_pbl_mix_humidity
    @test config.atmosphere_vertical_advection == :centered_full_temperature
    @test config.atmosphere_hyperdiffusion_hours == 4
    @test config.atmosphere_divergence_hyperdiffusion_hours == 4
    @test ExperimentConfig().atmosphere_divergence_hyperdiffusion_hours == 1
    spectral_grid = ReadyESM.SpeedyWeather.SpectralGrid(;
        NF = Float32, trunc = 7,
        Grid = ReadyESM.SpeedyWeather.RingGrids.FullGaussianGrid,
        nlayers = 8)
    horizontal_diffusion = ReadyESM._atmosphere_horizontal_diffusion(config, spectral_grid)
    @test Second(horizontal_diffusion.time_scale) == Second(14_400)
    @test Second(horizontal_diffusion.time_scale_div) == Second(14_400)
    @test config.atmosphere_large_scale_precipitation == :net_column_conservative
    @test config.atmosphere_saturation_law == :mixed_phase
    @test config.atmosphere_convection == :betts_miller_sigma_rh_v1
    @test config.atmosphere_convection_upper_relative_humidity == 0.5
    @test config.atmosphere_convection_lower_relative_humidity == 0.7
    @test config.atmosphere_convection_transition_top_sigma == 0.3
    @test config.atmosphere_convection_transition_bottom_sigma == 0.7
    @test config.atmosphere_convection_timescale_seconds == 7200
    @test config.atmosphere_initial_spectral_truncation
    @test !ExperimentConfig().atmosphere_initial_spectral_truncation
    @test config.forcing.radiation == :rrtmgp_all_sky
    @test config.forcing.solar_distance_correction
    @test !ExperimentConfig().forcing.solar_distance_correction
    @test config.forcing.co2_ppm == 420
    @test config.forcing.cloud_liquid_water_path_gm2 == 75
    @test config.forcing.cloud_ice_water_path_gm2 == 25
    @test config.forcing.aerosol_optical_depth_550nm == 0
end
println("PRODUCTION_CONFIG_TEST_PASS")
