module WriterGeometry

export geometry, define_geometry!, describe_variable!

"""Copy static native-grid coordinates. Cartesian ocean coordinates stay in metres."""
function geometry(x, y, z, area, active, longitude, latitude, weights; geographic=true)
    nx, ny = size(area)
    xx = ndims(x) == 1 ? repeat(reshape(x, :, 1), 1, ny) : x
    yy = ndims(y) == 1 ? repeat(reshape(y, 1, :), nx, 1) : y
    all(size(v) == (nx,ny) for v in (xx, yy, active)) || throw(DimensionMismatch("ocean geometry"))
    all(isfinite, xx) && all(isfinite, yy) && all(isfinite, z) || error("nonfinite ocean coordinates")
    all(isfinite, area) && all(area .>= 0) || error("invalid ocean cell areas")
    all(v -> v == 0 || v == 1, active) || error("ocean active mask is not binary")
    all(area[active .== 1] .> 0) || error("active ocean cell has zero area")
    geographic && !all(abs.(yy) .<= 90) && error("ocean latitude outside poles")
    length(longitude) == length(latitude) == length(weights) || throw(DimensionMismatch("atmosphere geometry"))
    all(isfinite, longitude) && all(isfinite, latitude) && all(abs.(latitude) .<= 90) || error("invalid atmosphere coordinates")
    all(isfinite, weights) && all(weights .> 0) || error("invalid atmosphere area weights")
    isapprox(sum(weights), 1; atol=2e-6, rtol=0) || error("atmosphere area weights must sum to one")
    isempty(z) && error("empty ocean vertical coordinate")
    all(diff(vec(z)) .> 0) || error("ocean z must increase upward")
    return (; x=Float64.(xx), y=Float64.(yy), z=Float64.(vec(z)),
        area=Float64.(area), active=Int8.(active), longitude=Float64.(vec(longitude)),
        latitude=Float64.(vec(latitude)), weights=Float64.(vec(weights)), geographic)
end

function define_geometry!(NC, ds, g)
    for (name,n) in (("ocean_i",size(g.area,1)), ("ocean_j",size(g.area,2)),
                     ("ocean_k",length(g.z)), ("atmosphere_point",length(g.weights)))
        NC.defDim(ds,name,n)
    end
    function variable(name,values,dims,units; attributes...)
        v=NC.defVar(ds,name,values,dims)
        v.attrib["units"]=units
        for (key,value) in attributes; v.attrib[String(key)]=value; end
    end
    xy=("ocean_i","ocean_j")
    if g.geographic
        variable("ocean_longitude",g.x,xy,"degrees_east"; standard_name="longitude")
        variable("ocean_latitude",g.y,xy,"degrees_north"; standard_name="latitude")
    else
        variable("ocean_x",g.x,xy,"m"; long_name="Cartesian ocean x coordinate")
        variable("ocean_y",g.y,xy,"m"; long_name="Cartesian ocean y coordinate")
    end
    variable("ocean_z",g.z,("ocean_k",),"m"; positive="up",
        long_name="reference ocean cell-centre vertical coordinate",
        comment="Reference grid coordinate, not the time-dependent physical height in z-star runs")
    variable("ocean_cell_area",g.area,xy,"m2"; long_name="native ocean surface-cell area")
    variable("ocean_active_mask",g.active,xy,"1";
        long_name="unique ocean surface cell for area reductions, 1 included and 0 excluded")
    variable("atmosphere_longitude",g.longitude,("atmosphere_point",),"degrees_east"; standard_name="longitude")
    variable("atmosphere_latitude",g.latitude,("atmosphere_point",),"degrees_north"; standard_name="latitude")
    variable("atmosphere_area_weight",g.weights,("atmosphere_point",),"1";
        long_name="normalized native Gaussian area weight per atmospheric point")
    ds.attrib["geometry_schema"]="native_writer_geometry_v1"
    ds.attrib["ocean_horizontal_coordinates"]=g.geographic ? "geographic" : "Cartesian"
    return nothing
end

const SPECS=Dict(
    "sst" => (:surface,"degree_Celsius","sea-surface temperature"),
    "sss" => (:surface,"psu","sea-surface salinity"),
    "ice_concentration" => (:surface,"1","sea-ice concentration"),
    "air_temperature_lowest" => (:atmosphere,"K","air temperature at the lowest model level"),
    "outgoing_longwave" => (:atmosphere,"W m-2","outgoing longwave radiation at the top of the atmosphere"),
    "rain_rate" => (:atmosphere,"m s-1","liquid precipitation rate, water equivalent"),
    "ocean_temperature" => (:volume,"degree_Celsius","ocean temperature"),
    "ocean_salinity" => (:volume,"psu","ocean salinity"),
    "global_sst" => (:scalar,"degree_Celsius","area-weighted global sea-surface temperature"),
    "global_outgoing_longwave" => (:scalar,"W m-2","area-weighted global outgoing longwave radiation"),
)

# Fields of kind :custom carry their own dimensions, in storage order, each with coordinate values and units, for
# example basin x level x latitude for an overturning streamfunction. A dimension is defined, with its coordinate
# variable, the first time a field uses it, and later fields must use the same coordinate values.
const CUSTOM_DIMS = Dict{String,Vector{Tuple{String,Vector{Float64},String}}}()
# Element type written for a field. Fields not listed keep Float64, so existing files are unchanged.
const STORAGE = Dict{String,DataType}()

function _custom_dims!(NC, ds, name, values)
    haskey(CUSTOM_DIMS, name) || error("missing dimensions for custom writer field: $name")
    spec = CUSTOM_DIMS[name]
    length(spec) == ndims(values) || throw(DimensionMismatch("custom field $name has $(ndims(values)) dimensions"))
    for (n, (d, coordinate, units)) in enumerate(spec)
        size(values, n) == length(coordinate) ||
            throw(DimensionMismatch("custom field $name: dimension $d has $(size(values, n)) values"))
        if haskey(ds.dim, d)
            Float64.(ds[d][:]) == coordinate || error("custom dimension $d already defined with other values")
        else
            NC.defDim(ds, d, length(coordinate))
            c = NC.defVar(ds, d, coordinate, (d,))
            c.attrib["units"] = units
        end
    end
    return Tuple(d for (d, _, _) in spec)
end

function describe_variable!(NC, ds, name, values, g; output_name=name, dims_override=nothing, mean=true)
    haskey(SPECS,name) || error("missing writer variable metadata: $name")
    kind,units,label=SPECS[name]
    dims = kind == :surface ? ("ocean_i","ocean_j") :
           kind == :volume ? ("ocean_i","ocean_j","ocean_k") :
           kind == :atmosphere ? ("atmosphere_point",) :
           kind == :custom && dims_override === nothing ? _custom_dims!(NC, ds, name, values) : ()
    if dims_override !== nothing
        kind == :scalar || error("dimension override requires a scalar diagnostic")
        dims=dims_override
    end
    expected=Tuple(ds.dim[d] for d in dims)
    size(values) == expected || throw(DimensionMismatch("writer field $name has $(size(values)), expected $expected"))
    v=NC.defVar(ds,output_name,get(STORAGE,name,Float64).(values),dims)
    v.attrib["units"]=units
    v.attrib["long_name"]=label
    v.attrib["sampling_statistic"]=mean ? "time mean" : "instantaneous launch state"
    if kind in (:surface,:volume)
        v.attrib["coordinates"]=g.geographic ? "ocean_longitude ocean_latitude" : "ocean_x ocean_y"
        kind == :volume && (v.attrib["coordinates"] *= " ocean_z")
        v.attrib["cell_measures"]="area: ocean_cell_area"
        v.attrib["ancillary_variables"]="ocean_active_mask"
    elseif kind == :atmosphere
        v.attrib["coordinates"]="atmosphere_longitude atmosphere_latitude"
        v.attrib["ancillary_variables"]="atmosphere_area_weight"
    end
    return v
end

end
