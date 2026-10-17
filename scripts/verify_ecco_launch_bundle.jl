using SHA

const CALENDAR_CENTRED_ECCO_BASE_FILES = Set((
    "CALENDAR_CENTRED_ECCO_1993_01_RECEIPT.json",
    "EVEL_1993_01.nc",
    "NVEL_1993_01.nc",
    "SALT_1993_01.nc",
    "SIarea_1993_01.nc",
    "SIheff_1993_01.nc",
    "SSH_1993_01.nc",
    "THETA_1993_01.nc",
))
const CALENDAR_CENTRED_ECCO_CACHE_FILES = Set((
    "EVEL_1993_01_inpainted.jld2",
    "NVEL_1993_01_inpainted.jld2",
    "SALT_1993_01_inpainted.jld2",
    "SSH_1993_01_inpainted.jld2",
    "THETA_1993_01_inpainted.jld2",
))
const CALENDAR_CENTRED_ECCO_FILES = union(
    CALENDAR_CENTRED_ECCO_BASE_FILES,
    CALENDAR_CENTRED_ECCO_CACHE_FILES,
)

function requires_calendar_centred_ecco_verification(
    config_name,
    directory,
    expected_directory,
)
    config_name == "production" && return true
    isdir(directory) && isdir(expected_directory) || return false
    return realpath(directory) == realpath(expected_directory)
end

function _calendar_centred_ecco_sha256(path)
    return open(path) do stream
        bytes2hex(sha256(stream))
    end
end

function _calendar_centred_ecco_manifest(path)
    isfile(path) || error("ECCO launch manifest is missing: $path")
    entries = Dict{String, String}()
    for (line_number, line) in enumerate(eachline(path))
        match_result = match(r"^([0-9a-f]{64})  ([A-Za-z0-9_.-]+)$", line)
        isnothing(match_result) && error(
            "malformed ECCO launch manifest line $line_number",
        )
        digest, filename = match_result.captures
        haskey(entries, filename) && error(
            "duplicate ECCO launch manifest entry: $filename",
        )
        entries[filename] = digest
    end
    Set(keys(entries)) == CALENDAR_CENTRED_ECCO_FILES || error(
        "ECCO launch manifest has an unexpected file set",
    )
    return entries
end

function verify_calendar_centred_ecco_bundle(
    directory,
    manifest_path;
    require_caches = false,
)
    isdir(directory) || error("ECCO launch directory is missing: $directory")
    entries = _calendar_centred_ecco_manifest(manifest_path)
    actual_files = Set(readdir(directory))
    valid_file_set = actual_files == CALENDAR_CENTRED_ECCO_BASE_FILES ||
        actual_files == CALENDAR_CENTRED_ECCO_FILES
    valid_file_set || error(
        "ECCO launch directory has an unexpected file set",
    )
    require_caches && actual_files != CALENDAR_CENTRED_ECCO_FILES && error(
        "ECCO launch directory lacks the complete qualified inpainting cache",
    )
    for filename in sort!(collect(actual_files))
        path = joinpath(directory, filename)
        islink(path) && error("ECCO launch file must not be a symbolic link: $path")
        actual = _calendar_centred_ecco_sha256(path)
        actual == entries[filename] || error(
            "ECCO launch SHA256 mismatch for $filename: expected " *
            "$(entries[filename]), found $actual",
        )
    end
    println(
        "CALENDAR_CENTRED_ECCO_VERIFICATION_PASS files=$(length(actual_files)) " *
        "caches=$(actual_files == CALENDAR_CENTRED_ECCO_FILES) " *
        "directory=$(abspath(directory))",
    )
    return nothing
end
