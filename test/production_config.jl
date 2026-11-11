using ReadyESM, Test, Dates, SHA

@testset "Released production configuration" begin
    config = load_config(joinpath(ReadyESM.PROJECT_ROOT, "config", "production.yml");
        check_input_files = false)
    @test config.device == :gpu
    @test config.duration_days == 730
    @test config.analysis_start_day == 365
    @test (config.truncation, config.nlayers) == (31, 27)
    @test config.ocean_grid == :tripolar_1degree
    @test (config.ocean_nlongitude, config.ocean_nlatitude, config.ocean_nlayers) == (360, 180, 60)
    @test config.ocean_dynamics
    @test config.sea_ice_dynamics
    @test config.sea_ice_lead_closing_thickness_m == 0.5
    @test config.land_model == :terrarium
    @test config.terrarium_precision == :float64
    @test config.terrarium_soil_layers == 16
    @test config.atmosphere_vertical_diffusion == :holtslag_boville_k_profile_v1
    @test config.atmosphere_pbl_diffusivity_scale == 1.0
    @test config.atmosphere_pbl_mix_momentum
    @test !config.atmosphere_pbl_mix_heat
    @test !config.atmosphere_pbl_mix_humidity
    @test config.atmosphere_dissipative_heating == :kinetic_energy_v1
    @test ExperimentConfig().atmosphere_dissipative_heating == :off
    @test config.atmosphere_vertical_advection == :centered_full_temperature
    @test config.atmosphere_hyperdiffusion_hours == 4
    @test config.atmosphere_divergence_hyperdiffusion_hours == 4
    @test ExperimentConfig().atmosphere_divergence_hyperdiffusion_hours == 1
    spectral_grid = ReadyESM.SpeedyWeather.SpectralGrid(;
        NF = Float32, trunc = 7,
        Grid = ReadyESM.SpeedyWeather.RingGrids.FullGaussianGrid,
        nlayers = 8)
    horizontal_diffusion = ReadyESM._atmosphere_horizontal_diffusion(config, spectral_grid)
    @test horizontal_diffusion isa ReadyESM.DissipativeHyperDiffusion
    @test horizontal_diffusion.heat
    plain_diffusion = ReadyESM._plain_horizontal_diffusion(horizontal_diffusion)
    @test Second(plain_diffusion.time_scale) == Second(14_400)
    @test Second(plain_diffusion.time_scale_div) == Second(14_400)
    @test config.atmosphere_large_scale_precipitation == :net_column_conservative
    @test !config.atmosphere_large_scale_melt_cooling
    @test !config.atmosphere_moist_mass_coupling
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
    @test config.forcing.well_mixed_gas_path == joinpath(ReadyESM.PROJECT_ROOT,
        "data", "forcing_1993", "wmghg_cmip6_ssp245_1850-2100.nc")
    @test config.forcing.well_mixed_gas_fixed_year == 1993
    @test config.forcing.ozone_path == joinpath(ReadyESM.PROJECT_ROOT,
        "data", "forcing_1993", "o3_cmip6_historical_199201-199512.nc")
    @test config.forcing.ozone_fixed_year == 0
    @test config.forcing.upper_radiative_cap_path == joinpath(ReadyESM.PROJECT_ROOT,
        "data", "upper_radiative_cap", "era5_1993_sampled_monthly_upper_cap.csv")
    @test bytes2hex(open(SHA.sha256, config.forcing.well_mixed_gas_path)) ==
        "d399d2d530ba946f4a53ec01a25caba956a8e2ac21646f5e02c2f45ef7f1b1b9"
    @test bytes2hex(open(SHA.sha256, config.forcing.ozone_path)) ==
        "fc687e57baba2df2e83182453827f08e1543131b5296602c298c72998b021d15"
    @test bytes2hex(open(SHA.sha256, config.forcing.upper_radiative_cap_path)) ==
        "5f07536991bd4a6a59536f7e5d86b58af803813c043d154d7a7d0d9184f5ca28"
    @test config.forcing.cloud_liquid_water_path_gm2 == 75
    @test config.forcing.cloud_ice_water_path_gm2 == 25
    @test config.forcing.cloud_relative_humidity_min == 0.3
    @test config.forcing.aerosol_optical_depth_550nm == 0
    @test config.forcing.aerosol_scheme == :cmip_prescribed
    @test config.forcing.solar_constant_entry == "picontrol"
    @test config.forcing.solar_constant_w_m2 == 1360.7470703125
    @test config.scored_outputs
    @test config.terrarium_surface_albedo == :era5_spectral
    @test config.terrarium_albedo_crossing_band == :equal_halves
    @test config.terrarium_albedo_cycle_sha256 ==
        "2f8e9193e2fea340e81dacc0e486c1d7de7f3b3a1b1e44dfc16ae3eef76f5c0c"
    @test config.terrarium_vegetation_cycle == :cover_lai
    @test config.terrarium_vegetation_cycle_sha256 ==
        "b80d7e01a4390340794030f1f924e5583d813fb163a1f9d875ae749091c13f60"
    @test ExperimentConfig().terrarium_surface_albedo == :constant
    @test ExperimentConfig().terrarium_vegetation_cycle == :off
    for (path, expected) in (
        (config.terrarium_albedo_cycle_path, config.terrarium_albedo_cycle_sha256),
        (config.terrarium_vegetation_cycle_path, config.terrarium_vegetation_cycle_sha256),
        (config.forcing.aerosol_cmip_stratospheric_path,
            "30e4e7383d6511c85c9aecdf84505bdb9f56d1aba717055aa0097c160daf607a"),
    )
        @test isfile(path)
        @test bytes2hex(open(sha256, path)) == expected
    end
end
println("PRODUCTION_CONFIG_TEST_PASS")
