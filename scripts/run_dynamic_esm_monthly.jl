#!/usr/bin/env julia

# Durable monthly scored diagnostics and a checkpoint paired with each sealed month.
# Resume starts in a fresh output directory and copies only the prefix named by
# the selected checkpoint receipt. An uncheckpointed later month stays untouched.

using ReadyESM, Dates, SHA, TOML
import CUDA

const R = ReadyESM
const O = R.Oceananigans
include(joinpath(@__DIR__, "monthly_windows", "monthly_windows.jl"))
const M = ReadyESMMonthlyWindows
include(joinpath(@__DIR__, "verify_ecco_launch_bundle.jl"))

sha(path) = bytes2hex(open(sha256, path))
const CHECKPOINT_SCHEMA = "readyesm_monthly_checkpoint_pair_v1"

function check_source()
    actual = R._dynamic_source_commit()
    expected = get(ENV, "READYESM_EXPECTED_MODEL_SOURCE", "")
    occursin(r"^[0-9a-f]{40}$", expected) || error("READYESM_EXPECTED_MODEL_SOURCE must be a 40-character commit")
    actual == expected || error("wrong ReadyESM source commit")
    isempty(readchomp(Cmd(["git", "-C", R.PROJECT_ROOT, "status",
                           "--porcelain", "--untracked-files=all"]))) ||
        error("ReadyESM source checkout is dirty")
    return actual
end

function check_production_ecco(config)
    expected = joinpath(R.PROJECT_ROOT, "data", "ecco_calendar_centred_1993_01")
    guard = requires_calendar_centred_ecco_verification(
        config.name, config.ecco_initial_conditions_directory, expected)
    if guard
        config.ocean_initial_conditions == :ecco_monthly_full_state ||
            error("production requires the qualified calendar-centred ocean launch")
        (config.ecco_year, config.ecco_month) == (1993, 1) ||
            error("production requires the January 1993 ECCO launch")
        realpath(config.ecco_initial_conditions_directory) == realpath(expected) ||
            error("production requires the qualified ECCO launch directory")
        verify_calendar_centred_ecco_bundle(
            config.ecco_initial_conditions_directory,
            joinpath(R.PROJECT_ROOT, "config", "ecco_calendar_centred_1993_01.sha256"))
    end
    return guard
end

function read_checkpoint_pair(path, config, source)
    isdir(path) || error("checkpoint directory missing: $path")
    record = TOML.parsefile(joinpath(path, "receipt.toml"))
    required = Set(("schema", "source_commit", "config_identity_sha256", "iteration",
                    "model_sha256", "month_receipts"))
    Set(keys(record)) == required || error("checkpoint receipt fields differ")
    record["schema"] == CHECKPOINT_SCHEMA || error("checkpoint receipt schema differs")
    record["source_commit"] == source || error("checkpoint source differs")
    record["config_identity_sha256"] == R._dynamic_config_identity_hash(config) ||
        error("checkpoint model configuration differs")
    iteration = record["iteration"]
    iteration isa Integer && iteration > 0 || error("invalid checkpoint iteration")
    checkpoint = joinpath(path, "model.jld2")
    sha(checkpoint) == record["model_sha256"] || error("checkpoint model digest differs")
    old_root = dirname(dirname(abspath(path)))
    dirname(abspath(path)) == joinpath(old_root, "checkpoints") ||
        error("checkpoint must be inside the old output's checkpoints directory")
    months = record["month_receipts"]
    months isa AbstractDict && !isempty(months) || error("empty month prefix")
    for (label, digest) in months
        occursin(r"^[0-9]{4}-[0-9]{2}$", label) || error("invalid month label: $label")
        receipt = joinpath(old_root, "months", label, "receipt.toml")
        sha(receipt) == digest || error("source month receipt digest differs: $label")
        month = TOML.parsefile(receipt)
        month["stop_iteration"] <= iteration || error("month later than checkpoint")
    end
    return (checkpoint=checkpoint, old_root=old_root, months=months, iteration=iteration)
end

function copy_prefix!(pair, output)
    target_root = joinpath(output, "months")
    mkpath(target_root)
    for label in sort!(collect(keys(pair.months)))
        source = joinpath(pair.old_root, "months", label)
        target = joinpath(target_root, label)
        !ispath(target) || error("month target exists: $target")
        cp(source, target)
        sha(joinpath(target, "receipt.toml")) == pair.months[label] ||
            error("copied month receipt differs: $label")
    end
end

function install_checkpoint_callback!(simulation, observer, config, source, prefix)
    root = joinpath(abspath(config.output_dir), "checkpoints")
    mkpath(root)
    written = Set{String}(prefix)
    identity = R._dynamic_config_identity_hash(config)
    function save_month_checkpoint(sim)
        iteration = Int(sim.model.clock.iteration)
        observer.callbacks == iteration ||
            error("month observer callback must precede checkpoint callback")
        fresh = [label for (label, record) in observer.records
                 if record["stop_iteration"] == iteration && !(label in written)]
        isempty(fresh) && return nothing
        length(fresh) == 1 || error("more than one month sealed at one step")
        label = only(fresh)
        target = joinpath(root, label)
        !ispath(target) || error("checkpoint target exists: $target")
        receipts = Dict{String, String}(
            key => String(value["receipt_sha256"]) for (key, value) in observer.records
            if value["stop_iteration"] <= iteration)
        Set(keys(receipts)) == union(written, Set((label,))) ||
            error("checkpoint month prefix differs")
        stage = mktempdir(root; prefix=".$label-incomplete-")
        model_path = joinpath(stage, "model.jld2")
        R.save_dynamic_restart_state(sim, model_path)
        record = Dict{String, Any}(
            "schema" => CHECKPOINT_SCHEMA,
            "source_commit" => source,
            "config_identity_sha256" => identity,
            "iteration" => iteration,
            "model_sha256" => sha(model_path),
            "month_receipts" => receipts,
        )
        open(joinpath(stage, "receipt.toml"), "w") do io
            TOML.print(io, record; sorted=true)
        end
        mv(stage, target)
        push!(written, label)
        println("V021_MONTH_CHECKPOINT_COMMITTED label=$label iteration=$iteration receipt_sha256=$(sha(joinpath(target, "receipt.toml")))")
        flush(stdout)
        return nothing
    end
    O.add_callback!(simulation, save_month_checkpoint, O.IterationInterval(1);
                    name=:readyesm_month_checkpoint_pair)
    steps_per_day = 86_400_000 ÷ observer.step_ms
    86_400_000 % observer.step_ms == 0 || error("monthly runner timestep must divide one day")
    function report_daily(sim)
        iteration = Int(sim.model.clock.iteration)
        iteration > 0 && iteration % steps_per_day == 0 || return nothing
        date = observer.run_start + Millisecond(iteration * observer.step_ms)
        println("V021_MONTHLY_PROGRESS date=$date iteration=$iteration")
        flush(stdout)
        return nothing
    end
    O.add_callback!(simulation, report_daily, O.IterationInterval(1);
                    name=:readyesm_monthly_progress)
    return nothing
end

function main(args=ARGS)
    length(args) in (2, 3) || error("usage: run_dynamic_esm_monthly.jl CONFIG.yml fresh|resume [OLD_CHECKPOINT_DIR]")
    config_path, mode = args[1:2]
    mode in ("fresh", "resume") || error("mode must be fresh or resume")
    length(args) == (mode == "fresh" ? 2 : 3) ||
        error("resume needs one old checkpoint directory")
    source = check_source()
    config = R.load_config(abspath(config_path))
    config.scored_outputs || error("monthly runner requires scored_outputs=true")
    output = abspath(config.output_dir)
    !ispath(output) || error("requires a fresh output directory: $output")
    if config.device == :gpu
        CUDA.functional() || error("GPU requested but CUDA is unavailable")
        CUDA.allowscalar(false)
        println("GPU: $(CUDA.name(CUDA.device()))")
    end
    guard = check_production_ecco(config)
    pair = mode == "resume" ? read_checkpoint_pair(abspath(args[3]), config, source) : nothing
    simulation = R.build_dynamic_esm(config)
    if guard
        verify_calendar_centred_ecco_bundle(
            config.ecco_initial_conditions_directory,
            joinpath(R.PROJECT_ROOT, "config", "ecco_calendar_centred_1993_01.sha256");
            require_caches=true)
    end
    if mode == "resume"
        copy_prefix!(pair, output)
    end
    observer = M.install_calendar_months!(simulation, config; resume=mode == "resume")
    if mode == "resume"
        Set(keys(observer.records)) == Set(keys(pair.months)) ||
            error("recovered month prefix differs from checkpoint receipt")
    end
    install_checkpoint_callback!(simulation, observer, config, source,
                                 mode == "resume" ? keys(pair.months) : ())
    if mode == "resume"
        M.restore_dynamic_restart_state!(simulation, pair.checkpoint)
        Int(simulation.model.clock.iteration) == pair.iteration ||
            error("restored model iteration differs from checkpoint receipt")
    end
    O.run!(simulation)
    inventory = M.save_windows!(observer, simulation, output)
    M.detach_completed!(observer, simulation)
    scored = R.save_scored_outputs(simulation, config)
    diagnostics = R.collect_dynamic_diagnostics(simulation, config)
    paths = R.save_dynamic_diagnostics(diagnostics, config; render_figure=false)
    R.validate_dynamic_diagnostics(diagnostics)
    R._write_dynamic_figure(paths.figure_path, diagnostics)
    println("V021_MONTHLY_RUN_PASS mode=$mode iteration=$(simulation.model.clock.iteration) months=$(inventory["committed_count"]) scored_hours=$(scored["monthly"]["state_samples"])")
    flush(stdout)
end

if abspath(PROGRAM_FILE) == (@__FILE__)
    main()
end
