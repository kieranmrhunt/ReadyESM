#!/usr/bin/env julia
# Check the committed dependency lock before Pkg.instantiate can rewrite it.
using TOML

root = dirname(@__DIR__)
project = TOML.parsefile(joinpath(root, "Project.toml"))
manifest = TOML.parsefile(joinpath(root, "Manifest.toml"))
name = project["name"]
entries = get(manifest["deps"], name, [])
length(entries) == 1 || error("Manifest must contain exactly one $name entry")
self = only(entries)
self["uuid"] == project["uuid"] || error("Manifest $name UUID differs from Project")
self["version"] == project["version"] || error("Manifest $name version differs from Project")
self["path"] == "." || error("Manifest $name must resolve to this source tree")
locked = self["deps"]
length(locked) == length(unique(locked)) || error("Manifest $name has duplicate dependencies")
expected = Set(keys(project["deps"]))
actual = Set(locked)
actual == expected || error("Manifest $name dependencies differ from Project: missing=$(sort!(collect(setdiff(expected, actual)))) extra=$(sort!(collect(setdiff(actual, expected))))")
for dep in expected
    haskey(manifest["deps"], dep) || error("Manifest lacks direct dependency $dep")
end
println("READY_ESM_MANIFEST_SELF_ENTRY_PASS version=$(project["version"]) deps=$(length(expected))")
