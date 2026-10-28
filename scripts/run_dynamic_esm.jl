#!/usr/bin/env julia

length(ARGS) == 1 || error("usage: run_dynamic_esm.jl CONFIG.yml")

using ReadyESM

import CUDA

config = load_config(only(ARGS))

include(joinpath(@__DIR__, "verify_ecco_launch_bundle.jl"))
expected_ecco_directory = joinpath(
    ReadyESM.PROJECT_ROOT,
    "data",
    "ecco_calendar_centred_1993_01",
)
guard_production_ecco = requires_calendar_centred_ecco_verification(
    config.name,
    config.ecco_initial_conditions_directory,
    expected_ecco_directory,
)
if guard_production_ecco
    config.ocean_initial_conditions == :ecco_monthly_full_state || error(
        "production requires the qualified calendar-centred ocean launch",
    )
    (config.ecco_year, config.ecco_month) == (1993, 1) || error(
        "production requires the qualified January 1993 ocean launch label",
    )
    realpath(config.ecco_initial_conditions_directory) ==
        realpath(expected_ecco_directory) || error(
            "production requires the qualified ECCO launch directory",
        )
    verify_calendar_centred_ecco_bundle(
        config.ecco_initial_conditions_directory,
        joinpath(
            ReadyESM.PROJECT_ROOT,
            "config",
            "ecco_calendar_centred_1993_01.sha256",
        ),
    )
end

if config.device == :gpu
    CUDA.functional() || error("device=gpu was requested but CUDA is not functional")
    CUDA.allowscalar(false)
    println("GPU: $(CUDA.name(CUDA.device()))")
end

println("Running $(config.name) with active $(config.forcing.co2_ppm) ppm CO2")
simulation = build_dynamic_esm(config)
if guard_production_ecco
    verify_calendar_centred_ecco_bundle(
        config.ecco_initial_conditions_directory,
        joinpath(
            ReadyESM.PROJECT_ROOT,
            "config",
            "ecco_calendar_centred_1993_01.sha256",
        );
        require_caches = true,
    )
end
ReadyESM.Oceananigans.run!(simulation)
if config.scored_outputs
    scored = ReadyESM.save_scored_outputs(simulation, config)
    println("SCORED_OUTPUTS_WRITTEN files=monthly.nc,wind.nc,toa_monthly.nc,eddy_6hourly.nc " *
            "hourly_samples=$(scored["monthly"]["state_samples"]) sixhourly_samples=$(scored["sixhourly"]["samples"])")
    flush(stdout)
end
diagnostics = collect_dynamic_diagnostics(simulation, config)
paths = save_dynamic_diagnostics(diagnostics, config; render_figure = false)
println("SIMULATION_COMPLETION_PASS iteration=$(simulation.model.clock.iteration) " *
        "diagnostics=$(paths.netcdf_path)")
flush(stdout)
validate_dynamic_diagnostics(diagnostics)
ReadyESM._write_dynamic_figure(paths.figure_path, diagnostics)
flush(stdout)
println("SCIENTIFIC_ACCEPTANCE_PASS")
println("Completed coupled iteration $(simulation.model.clock.iteration)")
println("Diagnostics: $(paths.netcdf_path)")
println("Figure: $(paths.figure_path)")
