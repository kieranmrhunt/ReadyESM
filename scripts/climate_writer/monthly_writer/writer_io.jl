# Bounded streaming output: NetCDF period files and writer-state checkpoints (Claude, prototype).
# The NCDatasets module is passed in by the caller (ReadyESM.NCDatasets), so this file has no package dependency.
module WriterIO

using Dates, SHA, Serialization
include(joinpath(@__DIR__, "writer_core.jl"))
using .WriterCore
const WC = WriterCore
include("writer_geometry.jl")
using .WriterGeometry

export WC, write_event, write_initial, save_state, load_state, file_name
export save_bound_state, load_bound_state, read_bound_snapshot, check_model_iteration

const CONTRACT = "sample = state after a coupled step at time t; period (start, end] containing t; the launch " *
                 "state is not in any mean"
const QUADRATURE = "right-endpoint rectangle rule, equal weights per coupled step"

file_name(e) = e.kind == :month ? "monthly_$(e.key[1])_$(lpad(e.key[2], 2, '0')).nc" :
               e.kind == :year ? "annual_$(e.key).nc" : error("unknown period kind $(e.kind)")

function fresh_paths(dir, name)
    path = joinpath(dir, name)
    tmp = path * ".tmp"
    ispath(path) && error("refusing to overwrite $path")
    ispath(tmp) && error("stale temporary file $tmp: inspect and remove it before resuming")
    return path, tmp
end

function common_attributes!(ds, cal, metadata)
    ds.attrib["writer_schema"] = "monthly_writer_prototype_v1"
    ds.attrib["calendar"] = "proleptic_gregorian"
    ds.attrib["origin"] = string(cal.origin)
    ds.attrib["step_seconds"] = cal.step_seconds
    ds.attrib["sampling_contract"] = CONTRACT
    ds.attrib["quadrature"] = QUADRATURE
    for (k, v) in sort(collect(metadata); by = first)
        ds.attrib[String(k)] = v
    end
    return ds
end

function define_arrays!(NC, ds, arrays; geometry=nothing, mean=true)
    for (name, value) in sort(collect(arrays); by = first)
        if geometry !== nothing
            describe_variable!(NC,ds,name,value,geometry; mean)
            continue
        end
        dims = ["$(name)_d$i" for i in 1:ndims(value)]
        for (d, n) in zip(dims, size(value))
            NC.defDim(ds, d, n)
        end
        NC.defVar(ds, name, Float64.(value), dims)
    end
    return ds
end

"""Write one completed period atomically. Files are never overwritten."""
function write_event(NC::Module, dir::AbstractString, e, cal::WC.Calendar, metadata::AbstractDict; geometry=nothing)
    mkpath(dir)
    path, tmp = fresh_paths(dir, file_name(e))
    NC.NCDataset(tmp, "c") do ds
        common_attributes!(ds, cal, metadata)
        ds.attrib["period"] = String(e.kind)
        ds.attrib["period_start"] = string(e.start)
        ds.attrib["period_end"] = string(e.stop)
        ds.attrib["first_sample_iteration"] = e.first_iteration
        ds.attrib["last_sample_iteration"] = e.last_iteration
        ds.attrib["first_sample_time"] = string(e.first_time)
        ds.attrib["last_sample_time"] = string(e.last_time)
        ds.attrib["sample_count"] = e.count
        ds.attrib["expected_sample_count"] = e.expected
        ds.attrib["complete"] = Int8(e.complete)
        ds.attrib["final_partial_flush"] = Int8(e.final)
        geometry === nothing || define_geometry!(NC,ds,geometry)
        define_arrays!(NC, ds, e.means; geometry)
        if e.kind == :month && !isempty(e.daily)
            NC.defDim(ds, "day", length(e.daily))
            NC.defVar(ds, "day_of_month", Int32[Dates.day(r.date) for r in e.daily], ("day",))
            for name in ("count", "expected", "first_iteration", "last_iteration")
                NC.defVar(ds, "day_$name", Int64[getproperty(r, Symbol(name)) for r in e.daily], ("day",))
            end
            for name in sort(collect(keys(first(e.daily).means)))
                values=Float64[r.means[name] for r in e.daily]
                if geometry === nothing
                    NC.defVar(ds,"daily_$name",values,("day",))
                else
                    describe_variable!(NC,ds,name,values,geometry;
                        output_name="daily_$name",dims_override=("day",))
                end
            end
        end
        ds.attrib["created_at_utc"] = string(now(UTC))
    end
    mv(tmp, path)
    return path
end

"""Write the launch state (iteration 0) separately from every period mean."""
function write_initial(NC::Module, dir::AbstractString, cal::WC.Calendar, fields::AbstractDict,
                       scalars::AbstractDict, metadata::AbstractDict; geometry=nothing)
    mkpath(dir)
    path, tmp = fresh_paths(dir, "initial.nc")
    NC.NCDataset(tmp, "c") do ds
        common_attributes!(ds, cal, metadata)
        ds.attrib["period"] = "initial"
        ds.attrib["sample_iteration"] = 0
        ds.attrib["sample_time"] = string(cal.origin)
        geometry === nothing || define_geometry!(NC,ds,geometry)
        define_arrays!(NC, ds, fields; geometry, mean=false)
        for (name, value) in sort(collect(scalars); by = first)
            values=fill(Float64(only(value)))
            if geometry === nothing
                NC.defVar(ds,"scalar_$name",values,())
            else
                describe_variable!(NC,ds,name,values,geometry; output_name="scalar_$name",mean=false)
            end
        end
        ds.attrib["created_at_utc"] = string(now(UTC))
    end
    mv(tmp, path)
    return path
end

"""Save the writer state next to a model checkpoint, with a SHA-256 sidecar. Never overwrites."""
function save_state(path::AbstractString, state::WC.State)
    ispath(path) && error("refusing to overwrite $path")
    bytes = let io = IOBuffer(); serialize(io, WC.snapshot(state)); take!(io) end
    write(path, bytes)
    write(path * ".sha256", bytes2hex(sha256(bytes)) * "\n")
    return path
end

"""Load and verify a writer state. The calendar must match and the model must be at the writer clock."""
function load_state(path::AbstractString, cal::WC.Calendar, model_iteration::Integer)
    bytes = read(path)
    bytes2hex(sha256(bytes)) == strip(read(path * ".sha256", String)) || error("writer state checksum differs: $path")
    return WC.restore_state(deserialize(IOBuffer(bytes)), cal, model_iteration)
end

include("writer_binding.jl")

end
