#!/usr/bin/env julia

using ReadyESM
using Dates

import Oceananigans

length(ARGS) in (1, 3, 4) || error(
    "usage: compare_dynamic_restart_boundaries.jl OUTPUT_DIRECTORY " *
    "[LEFT.jld2 RIGHT.jld2 [BASE_PATH]]",
)

directory = abspath(ARGS[1])
left_path = length(ARGS) >= 3 ? abspath(ARGS[2]) : joinpath(directory, "segment.jld2")
right_path = length(ARGS) >= 3 ? abspath(ARGS[3]) : joinpath(directory, "restored_boundary.jld2")
base_path = length(ARGS) == 4 ? ARGS[4] : "simulation"

const MAX_DIFFERENCES = 250

function array_difference_summary(left, right)
    # Materialize by iteration so arrays with offset axes compare by storage
    # order without asking Base.copyto! to reconcile incompatible indices.
    left_values = reshape([value for value in left], size(left))
    right_values = reshape([value for value in right], size(right))
    mask = .!isequal.(left_values, right_values)
    count = sum(mask)
    count == 0 && return nothing
    first_index = first(findall(mask))
    if eltype(left_values) <: Number && eltype(right_values) <: Number
        absolute_difference = abs.(left_values[mask] .- right_values[mask])
        finite = isfinite.(absolute_difference)
        max_abs = any(finite) ? maximum(absolute_difference[finite]) : NaN
        left_finite = isfinite.(left_values)
        right_finite = isfinite.(right_values)
        left_scale = any(left_finite) ? maximum(abs, left_values[left_finite]) : NaN
        right_scale = any(right_finite) ? maximum(abs, right_values[right_finite]) : NaN
        scale = max(left_scale, right_scale)
        relative_max = isfinite(max_abs) && isfinite(scale) ?
            max_abs / max(scale, floatmin(Float64)) : NaN
        return (; count, first_index, max_abs, scale, relative_max)
    end
    return (; count, first_index, max_abs = NaN, scale = NaN, relative_max = NaN)
end

function collect_differences!(differences, left, right, path)
    length(differences) >= MAX_DIFFERENCES && return differences
    if typeof(left) != typeof(right)
        push!(differences, "$path type $(typeof(left)) != $(typeof(right))")
        return differences
    elseif left isa Oceananigans.Fields.Field
        return collect_differences!(
            differences,
            parent(left),
            parent(right),
            "$path.data",
        )
    elseif left isa AbstractArray
        size(left) == size(right) || begin
            push!(differences, "$path size $(size(left)) != $(size(right))")
            return differences
        end
        summary = array_difference_summary(left, right)
        isnothing(summary) || push!(
            differences,
            "$path array count=$(summary.count)/$(length(left)) " *
            "first=$(summary.first_index) max_abs=$(summary.max_abs) " *
            "scale=$(summary.scale) relative_max=$(summary.relative_max)",
        )
        return differences
    elseif left isa NamedTuple
        keys(left) == keys(right) || begin
            push!(differences, "$path keys $(keys(left)) != $(keys(right))")
            return differences
        end
        for key in keys(left)
            key == :run_wall_time && continue
            collect_differences!(
                differences,
                getfield(left, key),
                getfield(right, key),
                "$path.$key",
            )
        end
        return differences
    elseif left isa ReadyESM.SpeedyWeather.Variables
        for name in propertynames(left)
            collect_differences!(
                differences,
                getfield(left, name),
                getfield(right, name),
                "$path.$name",
            )
        end
        return differences
    elseif left isa Base.RefValue
        return collect_differences!(differences, left[], right[], "$path[]")
    elseif left isa Number || left isa Symbol || left isa AbstractString ||
           left isa Dates.TimeType || left isa Dates.Period || isnothing(left)
        isequal(left, right) || push!(differences, "$path value $left != $right")
        return differences
    end

    names = fieldnames(typeof(left))
    isempty(names) && begin
        isequal(left, right) || push!(differences, "$path value differs")
        return differences
    end
    for name in names
        name == :run_wall_time && continue
        collect_differences!(
            differences,
            getfield(left, name),
            getfield(right, name),
            "$path.$name",
        )
    end
    return differences
end

left_loaded = Oceananigans.OutputWriters.load_checkpoint_state(
    left_path;
    base_path,
)
right_loaded = Oceananigans.OutputWriters.load_checkpoint_state(
    right_path;
    base_path,
)
left = base_path == "simulation" ? left_loaded.model : left_loaded
right = base_path == "simulation" ? right_loaded.model : right_loaded

differences = String[]
collect_differences!(differences, left, right, "model")
for (index, difference) in enumerate(differences)
    println("DYNAMIC_RESTART_BOUNDARY_DIFFERENCE index=$index $difference")
end
if isempty(differences)
    println("DYNAMIC_RESTART_BOUNDARY_EXACT_PASS")
else
    error(
        "restart boundaries differ in at least $(length(differences)) leaves" *
        (length(differences) == MAX_DIFFERENCES ? " (report capped)" : ""),
    )
end
