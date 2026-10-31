using NCDatasets, Dates, SHA
include("remap.jl")

"Read the validated 1993 ERA5 source, retaining four separate albedo fields."
function load_seasonal_source(path;expected_sha256)
    bytes2hex(open(sha256,path))==expected_sha256 || error("seasonal source hash differs")
    return NCDataset(path,"r") do ds
        time_name=haskey(ds,"valid_time") ? "valid_time" : "time"
        times=DateTime.(ds[time_name][:])
        times==[DateTime(1993,m,1) for m in 1:12] || error("wrong monthly source dates")
        latitude=Float64.(ds["latitude"][:]);longitude=Float64.(ds["longitude"][:])
        latitude==collect(90.0:-1.0:-90.0) || error("wrong source latitude")
        longitude==collect(0.0:359.0) || error("wrong source longitude")
        fields=Dict{Symbol,Array{Float64,3}}()
        for name in SeasonalLandRemapping.FIELD_NAMES
            variable=ds[String(name)]
            dimnames(variable)==("longitude","latitude",time_name) || error("wrong dimensions for $name")
            values=Float64.(coalesce.(variable[:,:,:],NaN))
            fields[name]=reverse(permutedims(values,(2,1,3));dims=1)
        end
        return SeasonalLandRemapping.source(reverse(latitude),longitude,fields;identity=expected_sha256)
    end
end
