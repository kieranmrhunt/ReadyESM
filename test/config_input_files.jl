@testset "Configuration inspection before input download" begin
    mktempdir() do directory
        source = joinpath(ReadyESM.PROJECT_ROOT, "config", "production.yml")
        settings = ReadyESM.YAML.load_file(source)
        for field in ("era5_pressure_levels_path", "era5_single_levels_path",
                      "era5_land_state_path")
            settings[field] = joinpath(directory, field * ".nc")
        end
        path = joinpath(directory, "experiment.yml")
        ReadyESM.YAML.write_file(path, settings)
        @test_throws ArgumentError load_config(path)
        inspected = load_config(path; check_input_files = false)
        @test inspected.atmosphere_initial_conditions == :era5_instantaneous
        @test inspected.terrarium_initial_conditions == :era5_instantaneous
        @test_throws ArgumentError ReadyESM.validate(inspected)
        # Construction rechecks inputs even when the configuration was inspected.
        @test_throws ArgumentError ReadyESM.build_dynamic_esm(inspected; steps = 1)

        # Inspection skips availability only; paths and physical settings are required.
        settings["era5_pressure_levels_path"] = ""
        ReadyESM.YAML.write_file(path, settings)
        @test_throws ArgumentError load_config(path; check_input_files = false)
        settings["era5_pressure_levels_path"] = joinpath(directory, "pressure.nc")
        settings["truncation"] = 1
        ReadyESM.YAML.write_file(path, settings)
        @test_throws ArgumentError load_config(path; check_input_files = false)
        settings["truncation"] = 31
        settings["ocean_initial_conditions"] = "ecco_monthly_full_state"
        settings["ecco_initial_conditions_directory"] = joinpath(directory, "ecco")
        ReadyESM.YAML.write_file(path, settings)
        inspected_ecco = load_config(path; check_input_files = false)
        @test_throws ArgumentError ReadyESM.validate(inspected_ecco)
    end
end
