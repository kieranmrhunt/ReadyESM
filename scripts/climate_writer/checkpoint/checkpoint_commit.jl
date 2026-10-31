module CheckpointCommit

using Dates, SHA, TOML

export at_month_end, checkpoint_path, commit!, read_commit, restore_prefix!

const PAYLOAD = ("model.jld2", "writer.jls", "writer.jls.sha256", "diagnostics.nc")
const OUTPUT_NAME = r"^(initial|monthly_[0-9]{4}_[0-9]{2}|annual_[0-9]{4})\.nc$"
digest(path) = bytes2hex(open(sha256, path))
sample_time(calendar, iteration) = calendar.origin + Second(iteration * calendar.step_seconds)
checkpoint_path(root, iteration) = joinpath(root, "checkpoints", "step_" * lpad(string(iteration), 10, '0'))

function at_month_end(calendar, iteration::Integer)
    iteration > 0 || return false
    time = sample_time(calendar, iteration)
    return time == DateTime(year(time), month(time), 1)
end

function fingerprint(path)
    isfile(path) && !islink(path) || error("expected a regular file: $path")
    filesize(path) > 0 || error("empty checkpoint file: $path")
    return Dict("bytes" => filesize(path), "sha256" => digest(path))
end

function check_file(path, expected)
    fingerprint(path) == expected || error("checkpoint file changed: $path")
    return nothing
end

function output_inventory(root)
    names = sort(filter(name -> occursin(OUTPUT_NAME, name), readdir(root)))
    "initial.nc" in names || error("missing launch output")
    return Dict(name => fingerprint(joinpath(root, name)) for name in names)
end

"""Publish one complete checkpoint directory by rename. Retain failed staging directories."""
function commit!(save_payload!, root, calendar, iteration::Integer, metadata::AbstractDict; reason)
    iteration > 0 || throw(ArgumentError("checkpoint iteration must be positive"))
    reason in ("month_end", "segment_end", "final") || throw(ArgumentError("unknown checkpoint reason"))
    reason == "month_end" && !at_month_end(calendar, iteration) && error("not a month boundary")
    target = checkpoint_path(root, iteration)
    ispath(target) && error("refusing to replace checkpoint: $target")
    prefix = output_inventory(root)
    parent = dirname(target)
    mkpath(parent)
    stage = mktempdir(parent; prefix = ".pending_$(iteration)_", cleanup = false)
    # The callback writes diagnostics, then a complete model/writer pair. It must
    # not mutate the writer's clock or publish output outside this staging path.
    save_payload!(stage)
    Set(readdir(stage)) == Set(PAYLOAD) || error("incomplete or unexpected checkpoint payload in $stage")
    files = Dict(name => fingerprint(joinpath(stage, name)) for name in PAYLOAD)
    for (name, expected) in prefix
        check_file(joinpath(root, name), expected)
    end
    record = Dict("schema" => "monthly_checkpoint_commit_v1", "iteration" => Int(iteration),
        "origin" => string(calendar.origin), "step_seconds" => calendar.step_seconds,
        "time" => string(sample_time(calendar, iteration)), "reason" => reason,
        "metadata" => Dict(String(k) => String(v) for (k, v) in metadata),
        "files" => files, "output_prefix" => prefix)
    receipt = joinpath(stage, "commit.toml")
    open(io -> TOML.print(io, record; sorted = true), receipt, "w")
    write(receipt * ".sha256", digest(receipt) * "\n")
    # Same-filesystem rename makes an incomplete pair unavailable under the
    # committed name. Checksums detect later damage. This is not an fsync claim.
    mv(stage, target)
    return target
end

function read_commit(path, calendar, metadata::AbstractDict)
    isdir(path) && !islink(path) || error("checkpoint directory is missing or linked")
    occursin(r"^step_[0-9]{10,}$", basename(path)) || error("checkpoint is not committed")
    receipt = joinpath(path, "commit.toml")
    isfile(receipt) && !islink(receipt) || error("commit receipt is missing")
    fingerprint(receipt * ".sha256")
    digest(receipt) == strip(read(receipt * ".sha256", String)) || error("commit receipt checksum differs")
    record = TOML.parsefile(receipt)
    record["schema"] == "monthly_checkpoint_commit_v1" || error("wrong commit schema")
    iteration = record["iteration"]
    iteration isa Int && iteration > 0 || error("invalid checkpoint iteration")
    basename(path) == basename(checkpoint_path("", iteration)) || error("checkpoint name/clock mismatch")
    record["origin"] == string(calendar.origin) && record["step_seconds"] == calendar.step_seconds ||
        error("checkpoint calendar differs")
    record["time"] == string(sample_time(calendar, iteration)) || error("checkpoint time differs")
    record["reason"] in ("month_end", "segment_end", "final") || error("unknown checkpoint reason")
    record["reason"] == "month_end" && !at_month_end(calendar, iteration) && error("not a month boundary")
    record["metadata"] == Dict(String(k) => String(v) for (k, v) in metadata) || error("checkpoint identity differs")
    Set(keys(record["files"])) == Set(PAYLOAD) || error("wrong checkpoint payload inventory")
    Set(readdir(path)) == Set((PAYLOAD..., "commit.toml", "commit.toml.sha256")) ||
        error("unexpected files in committed checkpoint")
    for (name, expected) in record["files"]
        check_file(joinpath(path, name), expected)
    end
    prefix = record["output_prefix"]
    haskey(prefix, "initial.nc") && all(name -> occursin(OUTPUT_NAME, name), keys(prefix)) ||
        error("invalid output prefix inventory")
    source_root = dirname(dirname(abspath(path)))
    for (name, expected) in prefix
        check_file(joinpath(source_root, name), expected)
    end
    return record
end

"""Copy only attested output into a fresh branch, never deleting later source output."""
function restore_prefix!(path, destination, calendar, metadata::AbstractDict)
    ispath(destination) && error("resume requires a fresh output directory")
    record = read_commit(path, calendar, metadata)
    source_root = dirname(dirname(abspath(path)))
    parent = dirname(abspath(destination))
    mkpath(parent)
    stage = mktempdir(parent; prefix = ".resume_pending_", cleanup = false)
    for (name, expected) in record["output_prefix"]
        source = joinpath(source_root, name)
        check_file(source, expected)
        cp(source, joinpath(stage, name))
        check_file(joinpath(stage, name), expected)
    end
    mv(stage, destination)
    return record
end

end
