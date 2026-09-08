module DarcyBoundaryTests

using ReadyESM, Test, KernelAbstractions
const R = ReadyESM
const T = R.Terrarium
const O = R.Oceananigans
const S = R.SpeedyWeather
const USE_GPU = get(ENV, "READYESM_DARCY_GPU", "false") == "true"
USE_GPU && R.CUDA.allowscalar(false)
device(x) = USE_GPU ? R.CUDA.CuArray(x) : x
sync() = USE_GPU ? R.CUDA.synchronize() : nothing

@kernel function sample_darcy!(selected, upstream, head, grid, psi, conductivity)
    i, k = @index(Global, NTuple)
    selected[i,k] = T.darcy_flux(i,1,k,grid,psi,conductivity)
    upstream[i,k] = invoke(T.darcy_flux, Tuple{Any,Any,Any,Any,Any,Any},
        i,1,k,grid,psi,conductivity)
    if k == grid.Nz + 1
        head[i,1] = psi[i,1,grid.Nz]
        head[i,2] = psi[i,1,grid.Nz+1]
    end
end

function sample(state, grid)
    fg = T.get_field_grid(grid)
    n = size(fg,1)
    NF = eltype(grid)
    selected, upstream = device(zeros(NF,n,fg.Nz+1)), device(zeros(NF,n,fg.Nz+1))
    head = device(zeros(NF,n,2))
    backend = get_backend(selected)
    sample_darcy!(backend)(selected,upstream,head,fg,
        state.pressure_head,state.hydraulic_conductivity; ndrange=size(selected))
    KernelAbstractions.synchronize(backend)
    return Array(selected),Array(upstream),Array(head)
end

function build_column(NF, geometry, model_kind; gradient=false)
    arch = USE_GPU ? T.GPU() : T.CPU()
    spacing = T.UniformSpacing(; Δz=0.1,N=10)
    grid = if geometry == :column
        T.ColumnGrid(arch,NF,spacing)
    else
        sg = S.SpectralGrid(; trunc=3,nlayers=2,Grid=R.RG.FullGaussianGrid)
        mask = ones(sg.grid) .> 0
        mask.data[1:2:end] .= false
        T.ColumnRingGrid(arch,NF,spacing,sg.grid,mask)
    end
    swrc = T.VanGenuchten(; α=NF(2),n=NF(2))
    hydraulics = T.ConstantSoilHydraulics(NF; swrc,
        unsat_hydraulic_cond=T.UnsatKVanGenuchten(NF))
    hydrology = T.SoilHydrology(NF,T.RichardsEq(); hydraulic_properties=hydraulics)
    soil = T.SoilEnergyWaterCarbon(NF; hydrology)
    model = if model_kind == :soil
        T.SoilModel(grid; soil,timestepper=T.ForwardEuler(NF))
    else
        T.LandModel(grid; soil,vegetation=nothing,snow=nothing,
            timestepper=T.ForwardEuler(NF))
    end
    # An explicit FluxBoundaryCondition has the same reflecting halo contract
    # as NoFlux. The gradient case is a separate allowed pressure-head BC.
    top = gradient ? O.GradientBoundaryCondition(NF(-0.25)) :
        O.FluxBoundaryCondition(NF(-1e-7))
    integrator = T.initialize(model; boundary_conditions=(; pressure_head=(; top)))
    state = integrator.state
    fg = T.get_field_grid(grid)
    n,Nz = size(fg,1),fg.Nz
    saturation = [NF(0.5)+NF(0.1)*NF((Nz-k+0.5)*0.1)
        for i in 1:n, j in 1:1, k in 1:Nz]
    T.set!(state.saturation_water_ice,saturation)
    T.set!(state.surface_excess_water,0)
    T.set!(state.temperature,15)
    if model_kind == :land
        for (name,value) in pairs((; air_temperature=15,air_pressure=101325,
            windspeed=2,specific_humidity=0.01,rainfall=1e-7,snowfall=0,
            surface_shortwave_down=0,surface_longwave_down=300))
            T.set!(getproperty(state.inputs,name),value)
        end
    end
    T.invclosure!(state,model)
    T.closure!(state,model)
    O.TimeSteppers.update_state!(integrator; compute_tendencies=true)
    sync()
    return integrator
end

@testset "Issue 1: actual plain/masked soil and land boundary lifecycle" begin
    for NF in (Float32,Float64), geometry in (:column,:masked_ring), kind in (:soil,:land)
        label = "$(NF)_$(geometry)_$(kind)"
        println("DARCY_CONSTRUCT $label device=$(USE_GPU ? "gpu" : "cpu")")
        flush(stdout)
        integrator = build_column(NF,geometry,kind)
        state, model = integrator.state,integrator.model
        fg = T.get_field_grid(model.grid)
        @test which(T.darcy_flux,Tuple{Int,Int,Int,typeof(fg),
            typeof(state.pressure_head),typeof(state.hydraulic_conductivity)}).module === T
        selected,upstream,head = sample(state,model.grid)
        @test head[:,1] == head[:,2]
        @test all(isfinite,upstream)
        @test all(iszero,upstream[:,end])
        @test selected == upstream
        @test any(!iszero,upstream[:,2:end-1])
        for step in 1:20
            T.timestep!(integrator,NF(60); finalize=false)
        end
        sync()
        O.TimeSteppers.update_state!(integrator; compute_tendencies=true)
        sync()
        selected,upstream,head = sample(state,model.grid)
        @test head[:,1] == head[:,2]
        @test all(iszero,upstream[:,end])
        @test selected == upstream
        prognostic = map(x -> Array(O.interior(x)),state.prognostic)
        @test all(all(isfinite,x) for x in values(prognostic))
        println("DARCY_REFLECTION_PASS $label native_steps=20 max_top_flux=$(maximum(abs,upstream[:,end]))")
        flush(stdout)
    end
end

@testset "Prescribed pressure gradient is not a no-flux boundary" begin
    for NF in (Float32,Float64), geometry in (:column,:masked_ring)
        integrator = build_column(NF,geometry,:soil; gradient=true)
        selected,upstream,head = sample(integrator.state,integrator.model.grid)
        @test all(head[:,2] .< head[:,1])
        @test all(upstream[:,end] .> 0)
        @test selected[:,1:end-1] == upstream[:,1:end-1]
        @test selected == upstream
    end
end

end

