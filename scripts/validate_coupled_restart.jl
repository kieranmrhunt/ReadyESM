using ReadyESM, Test
const R = ReadyESM
const O = R.Oceananigans
length(ARGS) == 2 || error("usage: validate_coupled_restart.jl CHECKPOINT NEW_OUTPUT")
checkpoint, output = abspath.(ARGS)
isfile(checkpoint) || error("missing full coupled checkpoint")
ispath(output) && error("refusing to overwrite restart evidence")
base = R.load_config(joinpath(R.PROJECT_ROOT, "config", "production.yml"))
settings = (; (name => getproperty(base, name) for name in fieldnames(typeof(base)))...)
config = R.ExperimentConfig(; merge(settings, (
    name = "public_native_runoff_fresh_restart",
    duration_days = 1 / 12,
    analysis_start_day = 0.0,
    output_dir = output,
))...)
R.CUDA.functional() || error("full coupled restart requires CUDA")
R.CUDA.allowscalar(false)
@test (config.truncation, config.nlayers) == (31, 27)
@test (config.ocean_nlongitude, config.ocean_nlatitude, config.ocean_nlayers) == (360, 180, 60)
println("PUBLIC_RUNOFF_FULL_RESTART_BUILD_BEGIN steps=8 restore_iteration=4")
flush(stdout)
simulation = R.build_dynamic_esm(config; steps=8)
@test Float64(simulation.Δt) == 900
println("PUBLIC_RUNOFF_FULL_RESTART_BUILD_PASS")
flush(stdout)

saved = O.OutputWriters.load_checkpoint_state(checkpoint; base_path="simulation")
@test saved.model.clock.iteration == 4
@test saved.model.clock.time == 3600
@test saved.model.land.schema == "readiesm_native_runoff_exchange_v1"
land = simulation.model.land
state = simulation.model.atmosphere.variables.prognostic.land.terrarium
@test !isnothing(land.runoff_exchange)
# Poison both owned stores before the normal full-model restore, so an omitted
# restore cannot accidentally pass because the freshly built fields are zero.
fill!(land.runoff_exchange.pending_depth, NaN)
fill!(land.surface_runoff, NaN)
R.restore_dynamic_restart_state!(simulation, checkpoint)
@test simulation.initialized
@test simulation.model.clock.iteration == 4
@test simulation.model.clock.time == 3600
@test O.prognostic_state(land) == saved.model.land
@test pointer(land.runoff_exchange.pending_depth) ==
    pointer(O.interior(state.runoff_pending_depth))
@test size(land.runoff_exchange.pending_depth) == size(O.interior(state.runoff_pending_depth))
@test strides(land.runoff_exchange.pending_depth) == strides(O.interior(state.runoff_pending_depth))
restored_path = joinpath(output, "restored_boundary.jld2")
R.save_dynamic_restart_state(simulation, restored_path)
println("PUBLIC_RUNOFF_FULL_RESTART_OWNED_STATE_PASS boundary=$restored_path")
flush(stdout)

# The saved boundary retains the original one-hour sampling horizon. Extend
# that control metadata only after saving the exact restored model state, so
# the atmosphere also records the endpoint of this additional hour.
atmosphere = simulation.model.atmosphere
atmosphere.variables.prognostic.clock.n_timesteps = 8
atmosphere.model.callbacks[:global_atmosphere_diagnostics].final_timestep = 8

# The restore helper has already reconciled and initialized the checkpoint.
# Native run! resets initialized=false, which would do that work a second
# time. Follow the existing initialized-restart verifier's stepping path.
simulation.stop_iteration = 8
simulation.running = true
while simulation.running
    O.time_step!(simulation)
end
for callback in values(simulation.callbacks)
    O.Simulations.finalize!(callback, simulation)
end
@test simulation.model.clock.iteration == 8
@test simulation.model.clock.time == 7200
for field in (land.runoff_exchange.pending_depth, land.surface_runoff)
    @test all(value -> isfinite(value) && value >= 0, Array(field))
end
@test pointer(land.runoff_exchange.pending_depth) ==
    pointer(O.interior(state.runoff_pending_depth))
diagnostics = R.collect_dynamic_diagnostics(simulation, config)
@test length(diagnostics.atmosphere_diagnostic_time_days) == 3
@test diagnostics.atmosphere_diagnostic_time_days[end] ≈ 1 / 12 rtol=eps(Float32)
paths = R.save_dynamic_diagnostics(diagnostics, config; render_figure=false)
R.save_dynamic_restart_state(simulation, joinpath(output, "hour2_restart.jld2"))
println("SIMULATION_COMPLETION_PASS iteration=8 diagnostics=$(paths.netcdf_path)")
flush(stdout)
R.validate_dynamic_diagnostics(diagnostics)
println("PUBLIC_RUNOFF_FULL_RESTART_EVOLUTION_PASS resumed_steps=4 full_climate=false trajectory_equality_not_tested=true")
flush(stdout)
