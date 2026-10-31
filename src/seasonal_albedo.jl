include("seasonal_albedo/adapter.jl")

"Land albedo supplied by the last spectral radiation solve."
struct HeldSpectralLandAlbedo{NF,RF} <: Terrarium.AbstractAlbedo{NF}
    initial_albedo::NF
    emissivity::NF
end
HeldSpectralLandAlbedo(::Type{NF},::Type{RF}) where {NF,RF} =
    HeldSpectralLandAlbedo{NF,RF}(NF(RF(COUPLED_LAND_SURFACE_ALBEDO)),NF(RF(COUPLED_SURFACE_EMISSIVITY)))
Adapt.adapt_structure(to,scheme::HeldSpectralLandAlbedo)=scheme
Terrarium.variables(scheme::HeldSpectralLandAlbedo)=(
    Terrarium.input(:albedo,Terrarium.XY();default=scheme.initial_albedo,
        bounds=Terrarium.UnitInterval,desc="Flux-equivalent land albedo held with radiation"),
)
@inline Terrarium.albedo(i,j,grid,fields,::HeldSpectralLandAlbedo)=@inbounds fields.albedo[i,j,1]
@inline Terrarium.emissivity(i,j,grid,fields,scheme::HeldSpectralLandAlbedo)=scheme.emissivity
Terrarium.compute_auxiliary!(state,grid,::HeldSpectralLandAlbedo,args...)=nothing

"Hold both surface components while a radiation solution is cached."
struct HeldSeasonalAtmosphereAlbedo{O} <: SpeedyWeather.AbstractAlbedo
    ocean_source::O
end
Adapt.@adapt_structure HeldSeasonalAtmosphereAlbedo
SpeedyWeather.initialize!(scheme::HeldSeasonalAtmosphereAlbedo,model::SpeedyWeather.PrimitiveEquation)=
    SpeedyWeather.initialize!(scheme.ocean_source,model)
@inline SpeedyWeather.albedo!(ij,field,vars,::HeldSeasonalAtmosphereAlbedo,model)=nothing

function _coupled_land_surface_energy_balance(::Type{FT},config::ExperimentConfig,::Type{RF}) where {FT,RF}
    config.terrarium_surface_albedo==:constant && return _coupled_land_surface_energy_balance(FT)
    return Terrarium.SurfaceEnergyBalance(FT;albedo=HeldSpectralLandAlbedo(FT,RF),
        skin_temperature=Terrarium.ImplicitSkinTemperature(FT;κₛ=FT(2)))
end
function _coupled_atmosphere_albedo(grid,config::ExperimentConfig)
    config.terrarium_surface_albedo==:constant && return _coupled_atmosphere_albedo(grid)
    return HeldSeasonalAtmosphereAlbedo(SpeedyWeather.OceanSeaIceAlbedo(grid))
end

struct SeasonalSpectralAlbedoController{M,I,NF} <: SpeedyWeather.AbstractCallback
    monthly::M
    indices::Vector{Int}
    weights::Vector{Float64}
    identity::I
    has_solved::Base.RefValue{Bool}
    direct::Matrix{NF}
    diffuse::Matrix{NF}
end

_has_seasonal_albedo(model)=hasproperty(model,:callbacks) && haskey(model.callbacks,:seasonal_albedo)
function _install_seasonal_albedo!(model,grid,config)
    config.terrarium_surface_albedo==:constant && return nothing
    model.albedo isa HeldSeasonalAtmosphereAlbedo || error("seasonal atmosphere boundary missing")
    land=model.land.model
    land.snow===nothing || throw(ArgumentError("spectral albedo currently requires snow-free land"))
    land.surface_energy_balance.albedo isa HeldSpectralLandAlbedo || error("held land boundary missing")
    source=load_seasonal_source(config.terrarium_albedo_cycle_path;expected_sha256=config.terrarium_albedo_cycle_sha256)
    longitude,latitude=SpeedyWeather.RingGrids.get_londlatds(grid.grid)
    indices=findall(vec(Bool.(Array(land.grid.mask.data))))
    fraction=vec(Array(model.land_sea_mask.mask.data))
    indices==findall(fraction.>0) || error("seasonal albedo must cover every positive land fraction")
    monthly=SeasonalLandRemapping.remap(source,longitude[indices],latitude[indices])
    solver=model.longwave_radiation.solver
    bounds=Array(RRTMGP.sw_band_bounds(solver))
    weights=SeasonalAlbedoAdapter.visible_weights(bounds;crossing=config.terrarium_albedo_crossing_band)
    buffer=IOBuffer()
    for values in (indices,Float64.(longitude[indices]),Float64.(latitude[indices]),
                   bounds,weights,(monthly.albedo[name] for name in SeasonalLandRemapping.ALBEDO_NAMES)...)
        write(buffer,reinterpret(UInt8,vec(values)))
    end
    mapping=bytes2hex(sha256(take!(buffer)))
    direct=RRTMGP.direct_sw_surface_albedo(solver)
    NF=eltype(direct)
    identity=(;schema=1,scheme=:era5_visible_nir_direct_diffuse,
        source_sha256=source.identity,mapping_sha256=mapping,
        calendar=:gregorian_month_centres,crossing_band=config.terrarium_albedo_crossing_band,
        exchange_precision=string(NF),land_precision=string(eltype(land.grid)))
    controller=SeasonalSpectralAlbedoController(monthly,indices,weights,identity,Ref(false),
        zeros(NF,size(direct)),zeros(NF,size(direct)))
    SpeedyWeather.add!(model.callbacks,:seasonal_albedo=>controller)
    @info "Seasonal spectral albedo mapped" columns=length(indices) fallback_count=count(monthly.fallback) maximum_fallback_degrees=maximum(monthly.distance) identity
    return controller
end

function _copy_held_albedo!(destination,values)
    length(destination)==length(values) || error("held albedo shape mismatch")
    copyto!(destination,reshape(eltype(destination).(values),size(destination)))
    KernelAbstractions.synchronize(KernelAbstractions.get_backend(destination))
    return nothing
end

KernelAbstractions.@kernel function _refresh_seasonal_ocean_albedo_kernel!(vars,scheme,model)
    ij=@index(Global,Linear)
    SpeedyWeather.albedo!(ij,vars.parameterizations.ocean.albedo,vars,scheme,model)
end
function _refresh_seasonal_ocean_albedo!(vars,model)
    backend=KernelAbstractions.get_backend(vars.parameterizations.ocean.albedo.data)
    _refresh_seasonal_ocean_albedo_kernel!(backend)(vars,model.albedo.ocean_source,model;ndrange=model.geometry.npoints)
    KernelAbstractions.synchronize(backend)
    return nothing
end

function SpeedyWeather.initialize!(controller::SeasonalSpectralAlbedoController,vars::SpeedyWeather.Variables,model::SpeedyWeather.AbstractModel)
    fill!(vars.parameterizations.land.albedo.data,COUPLED_LAND_SURFACE_ALBEDO)
    fill!(Oceananigans.interior(vars.prognostic.land.terrarium.inputs.albedo),COUPLED_LAND_SURFACE_ALBEDO)
    _refresh_seasonal_ocean_albedo!(vars,model)
    controller.has_solved[]=false
    return nothing
end
SpeedyWeather.callback!(::SeasonalSpectralAlbedoController,args...)=nothing
SpeedyWeather.finalize!(::SeasonalSpectralAlbedoController,args...)=nothing

function _prepare_seasonal_albedo!(radiation,vars,model)
    _has_seasonal_albedo(model) || return nothing
    controller=model.callbacks[:seasonal_albedo]
    _refresh_seasonal_ocean_albedo!(vars,model)
    fraction=vec(Array(model.land_sea_mask.mask.data))
    ocean=vec(Array(vars.parameterizations.ocean.albedo.data))
    b=SeasonalVegetation.bracket(vars.prognostic.clock.time)
    values=Dict{Symbol,Vector{Float64}}()
    for name in SeasonalLandRemapping.ALBEDO_NAMES
        monthly=controller.monthly.albedo[name]
        full=fill(NaN,length(fraction))
        full[controller.indices].=(1-b.weight).*monthly[:,b.lower].+b.weight.*monthly[:,b.upper]
        values[name]=full
    end
    SeasonalAlbedoAdapter.pack_albedos!(controller.direct,controller.diffuse,controller.weights,
        fraction,ocean,values[:aluvp],values[:alnip],values[:aluvd],values[:alnid])
    copyto!(RRTMGP.direct_sw_surface_albedo(radiation.solver),controller.direct)
    copyto!(RRTMGP.diffuse_sw_surface_albedo(radiation.solver),controller.diffuse)
    KernelAbstractions.synchronize(KernelAbstractions.get_backend(RRTMGP.direct_sw_surface_albedo(radiation.solver)))
    return nothing
end

function _finish_seasonal_albedo_solve!(radiation,vars,model)
    _has_seasonal_albedo(model) || return nothing
    controller=model.callbacks[:seasonal_albedo]
    fraction=vec(Array(model.land_sea_mask.mask.data))
    ocean=vec(Array(vars.parameterizations.ocean.albedo.data))
    held=vec(Array(vars.parameterizations.land.albedo.data))
    down=Array(RRTMGP.sw_flux_dn(radiation.solver))[1,:]
    up=Array(RRTMGP.sw_flux_up(radiation.solver))[1,:]
    NF=eltype(controller.direct)
    for col in eachindex(fraction)
        result=SeasonalAlbedoAdapter.equivalent_land_albedo(NF,fraction[col],ocean[col],down[col],up[col];previous=held[col])
        held[col]=result.albedo
    end
    _copy_held_albedo!(vars.parameterizations.land.albedo.data,held)
    _copy_held_albedo!(Oceananigans.interior(vars.prognostic.land.terrarium.inputs.albedo),held[controller.indices])
    mixture=(1 .- fraction).*ocean.+fraction.*held
    _copy_held_albedo!(vars.parameterizations.albedo.data,mixture)
    controller.has_solved[]=true
    return nothing
end

function _seasonal_albedo_restart_state(simulation)
    _has_seasonal_albedo(simulation.model) || return (;)
    controller=simulation.model.callbacks[:seasonal_albedo]
    solver=simulation.model.longwave_radiation.solver
    return (;seasonal_albedo=(;identity=controller.identity,has_solved=controller.has_solved[],
        direct=Array(RRTMGP.direct_sw_surface_albedo(solver)),
        diffuse=Array(RRTMGP.diffuse_sw_surface_albedo(solver))))
end

function _same_seasonal_albedo_identity(recorded,expected)
    recorded isa NamedTuple || return false
    length(recorded)==length(expected) || return false
    return all(name->hasproperty(recorded,name) &&
        isequal(getproperty(recorded,name),getproperty(expected,name)),keys(expected))
end

function _validate_seasonal_albedo_restart(simulation,state)
    enabled=_has_seasonal_albedo(simulation.model)
    present=hasproperty(state,:seasonal_albedo)
    enabled==present || throw(ArgumentError("checkpoint has different seasonal albedo mode"))
    enabled || return nothing
    controller=simulation.model.callbacks[:seasonal_albedo]
    record=state.seasonal_albedo
    required=(:identity,:has_solved,:direct,:diffuse)
    record isa NamedTuple && length(record)==length(required) &&
        all(name->hasproperty(record,name),required) ||
        throw(ArgumentError("invalid seasonal albedo checkpoint schema"))
    _same_seasonal_albedo_identity(record.identity,controller.identity) ||
        throw(ArgumentError("checkpoint has different seasonal albedo source, mapping or convention"))
    record.has_solved isa Bool || throw(ArgumentError("invalid seasonal albedo solve flag"))
    for name in (:direct,:diffuse)
        values=getproperty(record,name)
        size(values)==size(controller.direct) && eltype(values)==eltype(controller.direct) ||
            throw(ArgumentError("seasonal albedo checkpoint dimensions or precision differ"))
        all(SeasonalAlbedoAdapter.unit,values) || throw(ArgumentError("invalid checkpoint spectral albedos"))
    end
    return nothing
end

function _restore_seasonal_albedo_restart!(simulation,state)
    _has_seasonal_albedo(simulation.model) || return nothing
    _validate_seasonal_albedo_restart(simulation,state)
    record=state.seasonal_albedo
    solver=simulation.model.longwave_radiation.solver
    copyto!(RRTMGP.direct_sw_surface_albedo(solver),record.direct)
    copyto!(RRTMGP.diffuse_sw_surface_albedo(solver),record.diffuse)
    KernelAbstractions.synchronize(KernelAbstractions.get_backend(RRTMGP.direct_sw_surface_albedo(solver)))
    simulation.model.callbacks[:seasonal_albedo].has_solved[]=record.has_solved
    return nothing
end

function _seasonal_albedo_provenance(config)
    config.terrarium_surface_albedo==:constant && return (;)
    return (;terrarium_surface_albedo="era5_spectral",terrarium_albedo_source_sha256=config.terrarium_albedo_cycle_sha256,
        terrarium_albedo_calendar="gregorian_month_centres",terrarium_albedo_crossing_band=String(config.terrarium_albedo_crossing_band),
        terrarium_albedo_visible_threshold_cm1=14286.,terrarium_albedo_update="held_with_radiation_cache")
end
