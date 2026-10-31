
"Model-owned prescribed inputs. Buffers contain no independent physical memory."
struct SeasonalVegetationController{M,FT} <: SpeedyWeather.AbstractCallback
    monthly::M
    cover::Vector{FT}
    leaf_index::Vector{FT}
    mapping_sha256::String
end

function _install_seasonal_vegetation!(model,spectral_grid,config)
    config.terrarium_vegetation_cycle == :off && return nothing
    config.terrarium_vegetation_cycle == :cover_lai || error("unsupported vegetation cycle")
    source=load_seasonal_source(config.terrarium_vegetation_cycle_path;
        expected_sha256=config.terrarium_vegetation_cycle_sha256)
    longitude,latitude=SpeedyWeather.RingGrids.get_londlatds(spectral_grid.grid)
    mask=vec(Bool.(Array(model.land.model.grid.mask.data)))
    length(mask)==length(longitude)==length(latitude) || error("seasonal land grid mismatch")
    indices=findall(mask)
    monthly=SeasonalLandRemapping.remap(source,longitude[indices],latitude[indices])
    # Bind the exact native column order, coordinates and mapped leaf area.
    buffer=IOBuffer()
    for values in (indices,Float64.(longitude[indices]),Float64.(latitude[indices]),
                   monthly.cover,monthly.area_leaf_index)
        write(buffer,reinterpret(UInt8,vec(values)))
    end
    mapping=bytes2hex(sha256(take!(buffer)))
    FT=config.terrarium_precision == :float64 ? Float64 : Float32
    controller=SeasonalVegetationController(monthly,zeros(FT,length(indices)),zeros(FT,length(indices)),mapping)
    SpeedyWeather.add!(model.callbacks,:seasonal_vegetation=>controller)
    @info "Prescribed seasonal cover/LAI mapped" columns=length(indices) fallback_count=count(monthly.fallback) maximum_fallback_degrees=maximum(monthly.distance) source_sha256=source.identity mapping_sha256=mapping
    return controller
end

function _sample_seasonal_vegetation!(controller,date::DateTime)
    b=SeasonalVegetation.bracket(date)
    monthly=controller.monthly
    for i in eachindex(controller.cover)
        cover=(1-b.weight)*monthly.cover[i,b.lower]+b.weight*monthly.cover[i,b.upper]
        area=(1-b.weight)*monthly.area_leaf_index[i,b.lower]+b.weight*monthly.area_leaf_index[i,b.upper]
        controller.cover[i]=cover
        controller.leaf_index[i]=iszero(cover) ? 0 : area/cover
    end
    return nothing
end

function _update_seasonal_vegetation!(vars,model)
    haskey(model.callbacks,:seasonal_vegetation) || return nothing
    controller=model.callbacks[:seasonal_vegetation]
    _sample_seasonal_vegetation!(controller,vars.prognostic.clock.time)
    inputs=vars.prognostic.land.terrarium.inputs
    for (name,values) in ((:vegetation_fraction,controller.cover),(:leaf_area_index,controller.leaf_index))
        destination=Oceananigans.interior(getproperty(inputs,name))
        length(destination)==length(values) || error("seasonal vegetation column count differs")
        eltype(destination)==eltype(values) || error("seasonal vegetation precision differs")
        copyto!(destination,reshape(values,size(destination)))
    end
    # Existing root fractions, soil water and energy are untouched.
    return nothing
end

SpeedyWeather.initialize!(::SeasonalVegetationController,vars::SpeedyWeather.Variables,model::SpeedyWeather.AbstractModel) =
    _update_seasonal_vegetation!(vars,model)
SpeedyWeather.callback!(::SeasonalVegetationController,args...) = nothing
SpeedyWeather.finalize!(::SeasonalVegetationController,args...) = nothing

function _seasonal_vegetation_restart_state(simulation)
    callbacks=simulation.model.callbacks
    haskey(callbacks,:seasonal_vegetation) || return (;)
    controller=callbacks[:seasonal_vegetation]
    identity=(;schema=1,mode=:cover_lai,calendar=:gregorian_month_centres,
        source_sha256=controller.monthly.source_identity,mapping_sha256=controller.mapping_sha256)
    return (;seasonal_vegetation=identity)
end

function _validate_seasonal_vegetation_restart(simulation,state)
    expected=_seasonal_vegetation_restart_state(simulation)
    recorded=hasproperty(state,:seasonal_vegetation) ? (;seasonal_vegetation=state.seasonal_vegetation) : (;)
    isequal(expected,recorded) || throw(ArgumentError("checkpoint has different prescribed vegetation cycle or native mapping"))
    return nothing
end
