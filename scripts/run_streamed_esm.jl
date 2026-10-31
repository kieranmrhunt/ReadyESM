#!/usr/bin/env julia

# A read-only, opt-in monthly observer for long coupled runs. The ordinary
# run_dynamic_esm.jl path is unchanged. This runner is still experimental.

using ReadyESM, Dates, SHA, TOML
import CUDA

const R = ReadyESM
const O = R.Oceananigans
const SW = R.SpeedyWeather
const NC = R.NCDatasets
const CLIMATE_WRITER_DIR = joinpath(@__DIR__, "climate_writer")
const WRITER_LIBRARY = joinpath(CLIMATE_WRITER_DIR, "monthly_writer")
const WRITER_LIBRARY_SOURCE = "441547eacbabb0ccf99bb63c106ffb9f3b843eec"
include(joinpath(WRITER_LIBRARY, "writer_io.jl"))
include(joinpath(WRITER_LIBRARY, "writer_sampling.jl"))
include(joinpath(R.PROJECT_ROOT, "scripts", "verify_ecco_launch_bundle.jl"))
using .WriterIO
using .WriterSampling
const W = WriterIO.WC
const FIELD_BUNDLE = joinpath(CLIMATE_WRITER_DIR, "extended_fields")
include(joinpath(FIELD_BUNDLE, "sample.jl"))
const PROFILE_BUNDLE = joinpath(CLIMATE_WRITER_DIR, "profiles")
const ENERGY_BUNDLE = joinpath(CLIMATE_WRITER_DIR, "energy")
include(joinpath(PROFILE_BUNDLE, "sample.jl"))
include(joinpath(ENERGY_BUNDLE, "sample.jl"))
ClimateEnergy.register_metadata!(WriterIO.WriterGeometry.SPECS)
const CMIP_BUNDLE = joinpath(CLIMATE_WRITER_DIR, "cmip")
const CMIP_CORE_SHA256 = "ade96f1bd8a84d24d7b7b85a62530c82c3e72252fb0933b3aa85313e2db10ae6"
const CMIP_SAMPLER_SHA256 = "1ddaec0cc2efec7f717160de165f1c22938d3f823d14caf23e058a7df1ebe78a"
const CMIP_OCEAN_CORE_SHA256 = "6561ef7a547df4bce71248ab3c218c66397429581b52795b14db0419e1e78ed0"
const CMIP_OCEAN_SAMPLER_SHA256 = "94fe0b385693329e2feee0aeaf0c316f058665730e8876977ed073924f971d0d"
const CMIP_OVERTURNING_CORE_SHA256 = "2b9f7409d7341d77090e0f31debfdb45404478934011094799cadd4e5681915f"
const CMIP_OVERTURNING_SAMPLER_SHA256 = "292a5101d0e56129e8a1ed58a5ffc4dafd3856f539d6dd82135220b7a75feaf5"
bytes2hex(sha256(read(joinpath(CMIP_BUNDLE, "cmip_fields.jl")))) == CMIP_CORE_SHA256 || error("CMIP core changed")
bytes2hex(sha256(read(joinpath(CMIP_BUNDLE, "sample.jl")))) == CMIP_SAMPLER_SHA256 || error("CMIP sampler changed")
bytes2hex(sha256(read(joinpath(CMIP_BUNDLE, "ocean_fields.jl")))) == CMIP_OCEAN_CORE_SHA256 || error("CMIP ocean core changed")
bytes2hex(sha256(read(joinpath(CMIP_BUNDLE, "sample_ocean.jl")))) == CMIP_OCEAN_SAMPLER_SHA256 || error("CMIP ocean sampler changed")
bytes2hex(sha256(read(joinpath(CMIP_BUNDLE, "overturning.jl")))) == CMIP_OVERTURNING_CORE_SHA256 || error("CMIP overturning core changed")
bytes2hex(sha256(read(joinpath(CMIP_BUNDLE, "sample_overturning.jl")))) == CMIP_OVERTURNING_SAMPLER_SHA256 || error("CMIP overturning sampler changed")
include(joinpath(CMIP_BUNDLE, "sample.jl"))
include(joinpath(CMIP_BUNDLE, "sample_overturning.jl"))
include(joinpath(CMIP_BUNDLE, "sample_ocean.jl"))
const SNAPSHOT_HELPER = joinpath(CLIMATE_WRITER_DIR, "checkpoint", "snapshot.jl")
include(SNAPSHOT_HELPER)
include(joinpath(CLIMATE_WRITER_DIR, "checkpoint", "checkpoint_commit.jl"))
include(joinpath(CLIMATE_WRITER_DIR, "checkpoint", "guard.jl"))
using .ManualStepGuard
using .CheckpointCommit

function usage()
    error("usage: run_streamed_esm.jl CONFIG.yml off|full|segment|resume WRITER_DIR " *
          "[STOP_ITERATION] [COMMITTED_CHECKPOINT_DIRECTORY]")
end

function main(args=ARGS)
    3 <= length(args) <= 5 || usage()
    config_path, mode, output_directory = args[1:3]
    mode in ("off", "full", "segment", "resume") || usage()
    mode == "off" && length(args) != 3 && usage()
    mode == "full" && length(args) != 3 && usage()
    mode == "segment" && length(args) != 4 && usage()
    mode == "resume" && length(args) != 5 && usage()

    source = get(ENV, "READYESM_EXPECTED_MODEL_SOURCE", "")
    occursin(r"^[0-9a-f]{40}$", source) || error("expected model source is required")
    readchomp(`git -C $(R.PROJECT_ROOT) rev-parse HEAD`) == source || error("wrong model source")
    isempty(readchomp(`git -C $(R.PROJECT_ROOT) status --porcelain`)) || error("dirty model source")

    config_path = abspath(config_path)
    output_directory = abspath(output_directory)
    config = R.load_config(config_path)
    if config.device == :gpu
        CUDA.functional() || error("GPU requested but CUDA is unavailable")
        CUDA.allowscalar(false)
    end

    expected_ecco_directory = joinpath(R.PROJECT_ROOT, "data", "ecco_calendar_centred_1993_01")
    guard_production_ecco = requires_calendar_centred_ecco_verification(
        config.name, config.ecco_initial_conditions_directory, expected_ecco_directory,
    )
    if guard_production_ecco
        config.ocean_initial_conditions == :ecco_monthly_full_state ||
            error("production requires the qualified calendar-centred ECCO launch")
        (config.ecco_year, config.ecco_month) == (1993, 1) ||
            error("production requires the January 1993 ECCO launch")
        realpath(config.ecco_initial_conditions_directory) == realpath(expected_ecco_directory) ||
            error("production requires the qualified ECCO launch directory")
        verify_calendar_centred_ecco_bundle(
            config.ecco_initial_conditions_directory,
            joinpath(R.PROJECT_ROOT, "config", "ecco_calendar_centred_1993_01.sha256"),
        )
    end

    metadata = Dict(
        "source_commit" => source,
        "writer_library_source" => WRITER_LIBRARY_SOURCE,
        "extended_fields_sha256" => bytes2hex(sha256(read(joinpath(FIELD_BUNDLE, "fields.jl")))),
        "extended_sampler_sha256" => bytes2hex(sha256(read(joinpath(FIELD_BUNDLE, "sample.jl")))),
        "radiation_sampling" => "toa_* and surface radiation: latest postphysics cache; outgoing_longwave: completed-step applied flux",
        "launch_radiation" => "NaN if no solve has occurred, excluded from period means",
        "surface_radiation_domain" => "RRTMGP blended grid-box surface, not an individual ocean, ice or land tile",
        "ice_hemisphere_convention" => "latitude >= 0 is northern, latitude < 0 is southern",
        "configuration_sha256" => bytes2hex(sha256(read(config_path))),
        "writer_sha256" => bytes2hex(sha256(read(@__FILE__))),
        "step_guard_sha256" => bytes2hex(sha256(read(joinpath(CLIMATE_WRITER_DIR, "checkpoint", "guard.jl")))),
        "checkpoint_transaction_sha256" => bytes2hex(sha256(read(joinpath(CLIMATE_WRITER_DIR, "checkpoint", "checkpoint_commit.jl")))),
        "atmosphere_spatial_weights" => "model Gaussian ring quadrature divided by ring point count",
        "ocean_spatial_weights" => "unique active surface-cell area, duplicate center-fold storage excluded",
        "geometry_schema" => "native_writer_geometry_v1",
        "checkpoint_snapshot_sha256" => bytes2hex(sha256(read(SNAPSHOT_HELPER))),
        "checkpoint_snapshot_policy" => "output-only endpoint sample, persistent restart history unchanged",
        "energy_schema" => "native_latent_components_v1",
        "energy_core_sha256" => bytes2hex(sha256(read(joinpath(ENERGY_BUNDLE, "energy.jl")))),
        "energy_sampler_sha256" => bytes2hex(sha256(read(joinpath(ENERGY_BUNDLE, "sample.jl")))),
        "energy_ocean_timing" => "latest post-step interface cache, phase-weighted on unique wet ocean cells",
        "energy_land_timing" => "retained land substep-average humidity flux times its vaporisation constant, snow disabled",
        "energy_initial_availability" => "NaN, excluded from period means",
        "energy_aggregation" => "native component maps, no implicit land/ocean area renormalization",
        "profile_schema" => "climate_pressure_contributions_v1",
        "profile_core_sha256" => bytes2hex(sha256(read(joinpath(PROFILE_BUNDLE, "profiles.jl")))),
        "profile_sampler_sha256" => bytes2hex(sha256(read(joinpath(PROFILE_BUNDLE, "sample.jl")))),
        "profile_interpolation" => "linear in log pressure, no extrapolation",
        "profile_coverage" => "contribution and coverage averaged over every coupled step; score only valid_fraction == 1",
        "cmip_atmosphere_schema" => "cmip_atmosphere_samples_v1",
        "cmip_atmosphere_core_sha256" => CMIP_CORE_SHA256,
        "cmip_atmosphere_sampler_sha256" => CMIP_SAMPLER_SHA256,
        "cmip_atmosphere_timing" => "post-step; surface stress reconstructed from the applied u_prev/v_prev tendency inputs",
        "cmip_ocean_core_sha256" => CMIP_OCEAN_CORE_SHA256,
        "cmip_ocean_sampler_sha256" => CMIP_OCEAN_SAMPLER_SHA256,
        "cmip_overturning_core_sha256" => CMIP_OVERTURNING_CORE_SHA256,
        "cmip_overturning_sampler_sha256" => CMIP_OVERTURNING_SAMPLER_SHA256,
        "cmip_ocean_monthly_fields" => "zos, mlotst, 0-300, 300-700, 700-2000 m and below-2000 m heat content, msftmz, hfbasin and overturning face counts",
        "cmip_ocean_annual_fields" => "layer heat content, wet mask and native u and v velocities, alongside ocean temperature and salinity",
        "cmip_ocean_3d_cadence" => "annual only",
    )
    ispath(output_directory) && error("requires a fresh output directory: $output_directory")
    commit_directory = mode == "resume" ? abspath(args[5]) : nothing
    pending_writer = if mode == "resume"
        # Read the declared step only to form a calendar. All receipt/payload
        # hashes and the bound pair are verified before model construction.
        declared = TOML.parsefile(joinpath(commit_directory, "commit.toml"))
        provisional_calendar = W.Calendar(config.start_date, declared["step_seconds"])
        record = read_commit(commit_directory, provisional_calendar, metadata)
        payload = read_bound_snapshot(joinpath(commit_directory, "writer.jls"),
                                      joinpath(commit_directory, "model.jld2"), metadata)
        payload.model_iteration == record["iteration"] || error("commit/writer clock mismatch")
        payload
    else
        nothing
    end

    simulation = R.build_dynamic_esm(config)
    if guard_production_ecco
        verify_calendar_centred_ecco_bundle(
            config.ecco_initial_conditions_directory,
            joinpath(R.PROJECT_ROOT, "config", "ecco_calendar_centred_1993_01.sha256");
            require_caches = true,
        )
    end
    step_seconds = round(Int, Float64(simulation.Δt))
    Float64(simulation.Δt) == step_seconds || error("writer requires an integer-second coupled step")
    target_steps_float = config.duration_days * 86_400 / step_seconds
    target_steps = round(Int, target_steps_float)
    isapprox(target_steps_float, target_steps; atol = 1e-9, rtol = 0) ||
        error("requested duration is not an integer number of coupled steps")
    stop_iteration = mode in ("segment", "resume") ? parse(Int, args[4]) : target_steps
    0 < stop_iteration <= target_steps || error("stop iteration outside the configured run")
    mode == "segment" && stop_iteration == target_steps &&
        error("a segment must stop before the configured final iteration")

    earth = simulation.model
    ocean, ice, atmosphere = earth.ocean.model, earth.sea_ice.model, earth.atmosphere
    # Use the same normalized ring-area quadrature as the accepted radiation budget.
    # cos(latitude) per point is incorrect on the reduced octahedral Gaussian grid.
    atmosphere_weights = Float64.(Array(R._global_point_weights(atmosphere.model.spectral_grid)))
    ocean_active = unique_surface_mask(R._ocean_surface_active_mask(ocean.grid);
        right_center_folded=O.topology(ocean.grid,2)==O.RightCenterFolded)
    ocean_weights = Float64.(R._ocean_surface_cell_area(ocean.grid)) .* ocean_active
    geometry = if mode == "off"
        nothing
    else
        x,y,z = O.nodes(ocean.tracers.T)
        longitude,latitude = SW.RingGrids.get_londlatds(atmosphere.model.spectral_grid.grid)
        WriterIO.WriterGeometry.geometry(Array(x),Array(y),Array(z),
            R._ocean_surface_cell_area(ocean.grid),ocean_active,
            Array(longitude),Array(latitude),atmosphere_weights;
            # All three ReadyESM ocean constructors use spherical coordinates,
            # including the small "idealized" latitude-longitude grid.
            geographic=true)
    end

    function sample()
        temperature = Array(O.interior(ocean.tracers.T))
        salinity = Array(O.interior(ocean.tracers.S))
        surface_temperature = temperature[:, :, end]
        outgoing_longwave = vec(Array(atmosphere.variables.parameterizations.outgoing_longwave.data))
        monthly = Dict{String,Any}(
            "sst" => surface_temperature,
            "sss" => salinity[:, :, end],
            "ice_concentration" => Array(O.interior(ice.ice_concentration))[:, :, 1],
            "air_temperature_lowest" => Array(atmosphere.variables.grid.temperature.data)[:, end],
            "outgoing_longwave" => outgoing_longwave,
            "rain_rate" => vec(Array(atmosphere.variables.parameterizations.rain_rate.data)),
        )
        daily = Dict(
            "global_sst" => [weighted_mean(surface_temperature, ocean_weights)],
            "global_outgoing_longwave" => [weighted_mean(outgoing_longwave, atmosphere_weights)],
        )
        extra_monthly, extra_daily = additional_samples(R, earth, ocean_weights, atmosphere_weights; ocean_latitude=geometry.y)
        merge!(monthly, extra_monthly)
        merge!(monthly, climate_profile_samples(R, atmosphere))
        merge!(monthly, climate_energy_samples(R, earth, ocean_active))
        merge!(monthly, cmip_atmosphere_samples(R, earth))
        merge!(daily, extra_daily)
        annual = Dict("ocean_temperature" => temperature, "ocean_salinity" => salinity)
        cmip_ocean_monthly, cmip_ocean_annual = cmip_ocean_samples(
            R, earth; temperature=temperature, overturning=true)
        merge!(monthly, cmip_ocean_monthly)
        merge!(annual, cmip_ocean_annual)
        return monthly, daily, annual
    end

    calendar = W.Calendar(config.start_date, step_seconds)
    writer = nothing
    if mode == "resume"
        # Pair identity was checked before construction. Check the model's
        # actual calendar before restoration and its clock immediately after.
        writer = W.restore_state(pending_writer.state, calendar, pending_writer.model_iteration)
        R.restore_dynamic_restart_state!(simulation, joinpath(commit_directory, "model.jld2"))
        Float64(simulation.model.clock.time) == pending_writer.model_iteration * step_seconds ||
            error("restored model time differs from checkpoint iteration")
        check_model_iteration(writer, simulation.model.clock.iteration)
        simulation.model.clock.iteration < stop_iteration ||
            error("restored model is already at or past the requested stop iteration")
        restore_prefix!(commit_directory, output_directory, calendar, metadata)
        boundary = joinpath(output_directory, "restored_boundary_$(simulation.model.clock.iteration).jld2")
        ispath(boundary) && error("refusing to overwrite $boundary")
        R.save_dynamic_restart_state(simulation, boundary)
    else
        ispath(output_directory) && error("refusing to overwrite writer output: $output_directory")
        mkpath(output_directory)
        simulation.initialized = false
        simulation.running = true
        simulation.run_wall_time = 0.0
        O.Simulations.reset_nan_checker!(simulation)
        if mode != "off"
            writer = W.State(calendar)
            monthly, daily, _ = sample()
            write_initial(NC, output_directory, calendar, monthly, daily, metadata; geometry)
        end
    end

    simulation.running = true
    first_iteration = simulation.model.clock.iteration + 1
    sampling_seconds = Ref(0.0)
    checkpoint_seconds = Ref(0.0)
    function save_transaction(iteration, reason)
        check_model_iteration(writer, iteration)
        Float64(simulation.model.clock.time) == iteration * step_seconds ||
            error("model time differs from checkpoint iteration")
        path = commit!(output_directory, calendar, iteration, metadata; reason) do stage
            settings = (; (k => getproperty(config, k) for k in fieldnames(typeof(config)))...)
            snapshot_config = R.ExperimentConfig(; merge(settings, (; output_dir=stage))...)
            diagnostics = CheckpointSnapshot.collect_checkpoint_diagnostics(R, simulation, snapshot_config)
            R.save_dynamic_diagnostics(diagnostics, snapshot_config; render_figure=false)
            checkpoint = joinpath(stage, "model.jld2")
            R.save_dynamic_restart_state(simulation, checkpoint)
            save_bound_state(joinpath(stage, "writer.jls"), writer, checkpoint, metadata, iteration)
            # Keep all raw payload files in the failed stage if validation rejects them.
            R.validate_dynamic_diagnostics(diagnostics)
        end
        println("CHECKPOINT_COMMIT_PASS iteration=$iteration reason=$reason path=$path")
    end
    integration_seconds = @elapsed for iteration in first_iteration:stop_iteration
        advance_checked!(O.time_step!, s -> callback_nan_detected(O,s),
                         simulation, iteration, step_seconds, target_steps)
        if writer !== nothing
            sampling_seconds[] += @elapsed begin
                monthly, daily, annual = sample()
                for event in W.record!(writer, iteration, monthly, daily, annual)
                    write_event(NC, output_directory, event, calendar, metadata; geometry)
                end
            end
        end
        if writer !== nothing && (at_month_end(calendar, iteration) || iteration == stop_iteration)
            reason = at_month_end(calendar, iteration) ? "month_end" :
                     stop_iteration == target_steps ? "final" : "segment_end"
            checkpoint_seconds[] += @elapsed save_transaction(iteration, reason)
        end
        iteration % max(1, 86_400 ÷ step_seconds) == 0 &&
            println("STREAM_PROGRESS iteration=$iteration day=$(iteration * step_seconds / 86_400)")
    end

    # Save the resumable writer before final flushing. Partial final files are
    # deliberately outside its recorded prefix and are regenerated on a branch.
    if stop_iteration == target_steps && writer !== nothing
        for event in W.finalize!(writer)
            write_event(NC, output_directory, event, calendar, metadata; geometry)
        end
    end

    # Keep a full endpoint for the stronger coupled-state comparator, including
    # component velocities, lagged fluxes and callback histories.
    if stop_iteration == target_steps
        checkpoint = joinpath(output_directory, "model_endpoint_$(stop_iteration).jld2")
        ispath(checkpoint) && error("refusing to overwrite $checkpoint")
        R.save_dynamic_restart_state(simulation, checkpoint)
    end

    digest_array(values) = bytes2hex(sha256(reinterpret(UInt8, vec(Array(values)))))
    fingerprints = Dict(
        "iteration" => simulation.model.clock.iteration,
        "ocean_temperature" => digest_array(O.interior(ocean.tracers.T)),
        "ocean_salinity" => digest_array(O.interior(ocean.tracers.S)),
        "sea_ice_concentration" => digest_array(O.interior(ice.ice_concentration)),
        "atmosphere_temperature" => digest_array(atmosphere.variables.grid.temperature.data),
        "atmosphere_humidity" => digest_array(atmosphere.variables.grid.humidity.data),
    )
    fingerprint_path = joinpath(output_directory, "final_fingerprints_$(stop_iteration).toml")
    ispath(fingerprint_path) && error("refusing to overwrite $fingerprint_path")
    open(io -> TOML.print(io, fingerprints; sorted = true), fingerprint_path, "w")

    println("STREAM_RUN_PASS mode=$mode iteration=$(simulation.model.clock.iteration) " *
            "wall_seconds=$integration_seconds observer_seconds=$(sampling_seconds[]) checkpoint_seconds=$(checkpoint_seconds[])")
end

if abspath(PROGRAM_FILE) == (@__FILE__)
    main()
end
