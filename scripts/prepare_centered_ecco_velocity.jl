#!/usr/bin/env julia

# Prepare the two centered, CF-collocated ECCO4 velocity caches before a
# production model build. Run this on the eight-file calendar-centred input
# directory in isolation; the old staggered inpainting caches are not reused.
using ReadyESM, Dates, SHA

length(ARGS) == 1 || error(
    "usage: prepare_centered_ecco_velocity.jl CALENDAR_CENTRED_ECCO_DIRECTORY",
)
const R = ReadyESM
const NE = R.NumericalEarth
const directory = abspath(only(ARGS))
const manifest = joinpath(R.PROJECT_ROOT, "config", "ecco_calendar_centred_1993_01.sha256")
include(joinpath(@__DIR__, "verify_ecco_launch_bundle.jl"))

isdir(directory) || error("ECCO input directory is missing: $directory")
actual = Set(readdir(directory))
actual in (CALENDAR_CENTRED_ECCO_BASE_FILES,
           CALENDAR_CENTRED_ECCO_PREPARED_FILES,
           CALENDAR_CENTRED_ECCO_FILES) ||
    error("expected exact 8-file raw, 10-file prepared, or 13-file complete ECCO directory")
verify_calendar_centred_ecco_bundle(directory, manifest)
if actual == CALENDAR_CENTRED_ECCO_BASE_FILES
    expected = _calendar_centred_ecco_manifest(manifest)
    for (name, variable) in ((:u_velocity, "EVEL"), (:v_velocity, "NVEL"))
        metadata = NE.Metadatum(name; dataset=NE.ECCO4Monthly(),
            date=DateTime(1993, 1, 1), dir=directory)
        cache = R._ecco4_velocity_cache_path(metadata)
        ispath(cache) && error("refusing to replace existing velocity cache: $cache")
        field = R._ecco4_centered_velocity_field(metadata; require_cache=false)
        loaded = R._ecco4_centered_velocity_field(metadata)
        parent(field) == parent(loaded) || error("velocity cache readback differs: $cache")
        digest = bytes2hex(open(sha256, cache))
        digest == expected[basename(cache)] ||
            error("velocity cache SHA256 differs from qualified manifest: $cache")
        println("CENTERED_ECCO_VELOCITY_CACHE_PASS variable=$variable sha256=$digest")
        flush(stdout)
    end
end
verify_calendar_centred_ecco_bundle(directory, manifest)
println("CENTERED_ECCO_VELOCITY_PREPARATION_PASS files=$(length(readdir(directory)))")
