using ReadyESM, Test
const R = ReadyESM
const O = R.Oceananigans
length(ARGS) == 1 || error("usage: validate_coupled_smoke.jl OUTPUT")
output = abspath(only(ARGS))
ispath(output) && error("refusing to overwrite smoke evidence")
base = R.load_config(joinpath(R.PROJECT_ROOT, "config", "production.yml"))
settings = (; (name => getproperty(base, name) for name in fieldnames(typeof(base)))...)
config = R.ExperimentConfig(; merge(settings, (
    name = "public_native_runoff_four_step",
    duration_days = 1 / 24,
    analysis_start_day = 0.0,
    output_dir = output,
))...)
R.CUDA.functional() || error("coupled smoke requires CUDA")
R.CUDA.allowscalar(false)
@test (config.truncation, config.nlayers) == (31, 27)
@test (config.ocean_nlongitude, config.ocean_nlatitude, config.ocean_nlayers) == (360, 180, 60)
println("PUBLIC_RUNOFF_COUPLED_BUILD_BEGIN steps=4 full_climate=false")
flush(stdout)
simulation = R.build_dynamic_esm(config; steps=4)
@test Float64(simulation.Δt) == 900
@test !isnothing(simulation.model.land.runoff_exchange)
println("PUBLIC_RUNOFF_COUPLED_BUILD_PASS")
flush(stdout)
O.run!(simulation)
@test simulation.model.clock.iteration == 4
@test simulation.model.clock.time == 3600
diagnostics = R.collect_dynamic_diagnostics(simulation, config)
paths = R.save_dynamic_diagnostics(diagnostics, config; render_figure=false)
checkpoint = joinpath(output, "hour1_restart.jld2")
ispath(checkpoint) && error("refusing to overwrite $checkpoint")
R.save_dynamic_restart_state(simulation, checkpoint)
@test O.prognostic_state(simulation.model.land).schema == "readiesm_native_runoff_exchange_v1"
println("SIMULATION_COMPLETION_PASS iteration=4 diagnostics=$(paths.netcdf_path) checkpoint=$checkpoint")
flush(stdout)
R.validate_dynamic_diagnostics(diagnostics)
println("PUBLIC_RUNOFF_COUPLED_SMOKE_PASS seconds=3600 full_climate=false")
flush(stdout)
