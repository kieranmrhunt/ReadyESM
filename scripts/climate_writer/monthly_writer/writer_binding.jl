# Included inside WriterIO. A bound state belongs to one exact model checkpoint
# and one source/configuration/writer identity. Paths may change, bytes may not.

function _binding_metadata(metadata::AbstractDict)
    result = Dict(String(k) => String(v) for (k, v) in metadata)
    for key in ("source_commit", "configuration_sha256", "writer_sha256")
        haskey(result, key) && !isempty(result[key]) || error("missing writer binding metadata: $key")
    end
    return result
end

function check_model_iteration(state::WC.State, model_iteration::Integer)
    state.last_iteration == model_iteration ||
        error("writer clock $(state.last_iteration) differs from model iteration $model_iteration")
    return state
end

"""Save a writer state bound to exact model-checkpoint bytes and run identity. Never overwrites."""
function save_bound_state(path::AbstractString, state::WC.State, checkpoint::AbstractString,
                          metadata::AbstractDict, model_iteration::Integer)
    check_model_iteration(state, model_iteration)
    identity = _binding_metadata(metadata)
    isfile(checkpoint) || error("model checkpoint is missing: $checkpoint")
    for target in (path, path * ".sha256", path * ".tmp", path * ".sha256.tmp")
        ispath(target) && error("refusing to overwrite writer checkpoint artifact: $target")
    end
    payload = (
        schema = "monthly_writer_bound_state_v1",
        model_checkpoint_sha256 = bytes2hex(open(sha256, checkpoint)),
        model_iteration = Int(model_iteration),
        metadata = identity,
        state = WC.snapshot(state),
    )
    bytes = let io = IOBuffer()
        serialize(io, payload)
        take!(io)
    end
    # Incomplete pairs are rejected on load. Leave temporary files for inspection
    # if a write fails rather than silently accepting or replacing a partial state.
    write(path * ".tmp", bytes)
    write(path * ".sha256.tmp", bytes2hex(sha256(bytes)) * "\n")
    mv(path * ".tmp", path)
    mv(path * ".sha256.tmp", path * ".sha256")
    return path
end

"""Verify the saved bytes and pair identity before constructing or restoring the model."""
function read_bound_snapshot(path::AbstractString, checkpoint::AbstractString, metadata::AbstractDict)
    identity = _binding_metadata(metadata)
    bytes = read(path)
    bytes2hex(sha256(bytes)) == strip(read(path * ".sha256", String)) ||
        error("writer state checksum differs: $path")
    payload = deserialize(IOBuffer(bytes))
    hasproperty(payload, :schema) && payload.schema == "monthly_writer_bound_state_v1" ||
        error("writer state has no model-checkpoint binding: $path")
    payload.metadata == identity || error("writer source, configuration or sampling identity differs")
    isfile(checkpoint) || error("model checkpoint is missing: $checkpoint")
    bytes2hex(open(sha256, checkpoint)) == payload.model_checkpoint_sha256 ||
        error("writer state belongs to a different model checkpoint")
    return payload
end

"""Verify a writer/model checkpoint pair and calendar before restoring the model.
After model restoration, call `check_model_iteration` with the actual model clock."""
function load_bound_state(path::AbstractString, checkpoint::AbstractString,
                          metadata::AbstractDict, calendar::WC.Calendar)
    payload = read_bound_snapshot(path, checkpoint, metadata)
    return WC.restore_state(payload.state, calendar, payload.model_iteration)
end
