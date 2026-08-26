#!/usr/bin/env julia

using ReadyESM

import CUDA

config_path = isempty(ARGS) ?
    joinpath(ReadyESM.PROJECT_ROOT, "config", "dynamic_esm_smoke.yml") : ARGS[1]
config = load_config(config_path)

if config.device == :gpu
    CUDA.functional() || error("device=gpu was requested but CUDA is not functional")
    CUDA.allowscalar(false)
    println("GPU: $(CUDA.name(CUDA.device()))")
end

println("Running $(config.name) with active $(config.forcing.co2_ppm) ppm CO2")
simulation = run_dynamic_esm!(config)
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
