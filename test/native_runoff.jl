module NativeRunoffTests
using ReadyESM, Test, Dates, JLD2
const R = ReadyESM
const O = R.Oceananigans
const T = R.Terrarium
const USE_GPU = get(ENV,"READYESM_RUNOFF_GPU","false") == "true"
USE_GPU && R.CUDA.allowscalar(false)
device(x) = USE_GPU ? R.CUDA.CuArray(x) : x
sync() = USE_GPU ? R.CUDA.synchronize() : nothing
host(field) = Float64.(Array(O.interior(field)))

function synthetic_router(rate, pending)
    return R.TerrariumRunoffLand(rate,Int[],Float32[],Int[],Int[],[1],
        length(rate),0,16,:equal_area_flux,1.,nothing,(;pending_depth=pending))
end

@testset "Owned native-runoff amounts and restart" begin
    for FT in (Float32,Float64)
        amounts = reshape(FT[0,1e-9,0.004455836,0.02],:,1,1)
        pending = device(copy(amounts))
        rate = device(zero.(amounts))
        land = synthetic_router(rate,pending)
        O.TimeSteppers.time_step!(land,900)
        @test all(Array(rate) .>= 0)
        @test all(Array(pending) .>= 0)
        @test Float64.(Array(rate)).*900 .+ Float64.(Array(pending)) ≈
            Float64.(amounts) rtol=8eps(FT) atol=0
        saved = O.prognostic_state(land)
        @test saved.schema == "readiesm_native_runoff_exchange_v1"
        @test_throws ErrorException O.restore_prognostic_state!(land,nothing)
        @test_throws ErrorException O.restore_prognostic_state!(land,(;schema="unknown"))
        mktempdir() do directory
            path = joinpath(directory,"runoff.jld2")
            jldsave(path;state=saved)
            fill!(pending,NaN)
            fill!(rate,NaN)
            O.restore_prognostic_state!(land,load(path,"state"))
            @test O.prognostic_state(land) == saved
        end
        @test_throws ArgumentError O.TimeSteppers.time_step!(land,0)
        @test_throws ArgumentError O.TimeSteppers.time_step!(land,NaN)
        copyto!(pending,fill(FT(-1),size(amounts)))
        O.TimeSteppers.time_step!(land,900)
        @test all(isnan,Array(rate))
        @test all(==(-1),Array(pending))
        O.Simulations.reset_clock!(land)
        sync()
        @test all(iszero,Array(pending))
        @test all(iszero,Array(rate))
    end
    legacy = R.TerrariumRunoffLand(zeros(1,1,1),Int[],Float32[],Int[],Int[],[1],
        1,0,16,:equal_area_flux,1.)
    @test isnothing(O.prognostic_state(legacy))
    @test O.restore_prognostic_state!(legacy,nothing) === legacy
    @test isnothing(O.TimeSteppers.time_step!(legacy,900))
end

@testset "Production land pulse reaches the actual ocean router" begin
    config = R.ExperimentConfig(device=USE_GPU ? :gpu : :cpu,truncation=7,nlayers=8,
        land_model=:terrarium,terrarium_soil_layers=3,terrarium_max_layer_thickness_m=0.5,
        terrarium_timestep_seconds=60,forcing=R.ForcingConfig(radiation=:rrtmgp_clear_sky))
    simulation = R.build_rrtmgp(config)
    atmosphere = simulation.model
    land_model = atmosphere.land
    state = simulation.variables.prognostic.land.terrarium
    model = land_model.model
    FT = eltype(state)
    @test hasproperty(state,:runoff_pending_depth)
    for (name,value) in pairs((;air_temperature=10,specific_humidity=0.0075,
        air_pressure=101325,windspeed=1,rainfall=0,snowfall=0,
        surface_shortwave_down=0,surface_longwave_down=350))
        T.set!(getproperty(state.inputs,name),value)
    end
    T.set!(state.temperature,10)
    T.invclosure!(state,model)
    T.set!(state.saturation_water_ice,1)
    T.set!(state.surface_excess_water,0.02)
    T.closure!(state,model)
    T.compute_auxiliary!(state,model)
    saved = deepcopy(state)
    integrator = T.ModelIntegrator(state.clock,model,T.InputSources(FT),state,land_model.initializers)
    native = zeros(size(host(state.surface_runoff)))
    for half in 1:2
        for step in 1:8
            dt = step == 8 ? 30. : 60.
            T.timestep!(integrator,dt;finalize=false)
            native .+= dt .* host(state.surface_runoff)
        end
        T.compute_auxiliary!(state,model)
    end
    physical = map(field->Array(O.interior(field)),state.prognostic)
    endpoint = 900 .* host(state.surface_runoff)
    copyto!(state,saved)
    for half in 1:2
        R._run_terrarium_land!(integrator,Second(450),60.)
    end
    for name in keys(state.prognostic)
        @test Array(O.interior(getproperty(state.prognostic,name))) == getproperty(physical,name)
    end
    @test maximum(native) > 1e-4
    @test maximum(abs.(endpoint.-native)) > 0.01maximum(native)
    pending = O.interior(state.runoff_pending_depth)
    @test host(state.runoff_pending_depth) ≈ native rtol=16eps(FT) atol=0
    grid_cpu = O.TripolarGrid(O.CPU(),Float32;size=(16,12,2),halo=(3,3,3),z=(-100,0))
    arch = USE_GPU ? O.GPU() : O.CPU()
    grid = O.on_architecture(arch,grid_cpu)
    router = R._build_terrarium_runoff_land(simulation,(;model=(;grid)))
    @test size(router.runoff_exchange.pending_depth) == size(pending)
    @test strides(router.runoff_exchange.pending_depth) == strides(pending)
    @test pointer(router.runoff_exchange.pending_depth) == pointer(pending)
    saved_native = deepcopy(state)
    boundary = O.prognostic_state(router)
    O.TimeSteppers.time_step!(router,900)
    @test Float64.(Array(router.surface_runoff)).*900 .+ Array(pending) ≈
        native rtol=16eps(FT) atol=0
    expected = O.prognostic_state(router)
    flux = O.Field{O.Center,O.Center,Nothing}(grid)
    R._scatter_terrarium_runoff_state!(flux,router)
    O.Architectures.synchronize(arch)
    @test all(isfinite,Array(O.interior(flux)))
    weights = Float64.(Array(R._global_point_weights(atmosphere.spectral_grid)))
    fraction = Float64.(Array(atmosphere.land_sea_mask.mask.data))
    points = findall(fraction .> 0)
    area = weights[points].*fraction[points].*(4pi*atmosphere.planet.radius^2)
    exported = R._budget_scalar(R._budget_integral(flux))*900
    @test exported + 1000sum(area.*vec(Array(pending))) ≈
        1000sum(area.*vec(native)) rtol=2e-6
    before_fill = Array(O.interior(flux))
    O.fill_halo_regions!(flux)
    O.Architectures.synchronize(arch)
    @test Array(O.interior(flux)) == before_fill
    snapshot = R._runoff_surface_grid_matrix(flux)
    @test all(iszero,snapshot[9:16,12])
    @test snapshot[1:8,12] == Float64.(before_fill[1:8,12,1])
    @test Array(O.interior(flux)) == before_fill
    fill!(pending,NaN)
    fill!(router.surface_runoff,NaN)
    copyto!(state,saved_native)
    O.restore_prognostic_state!(router,boundary)
    @test pointer(router.runoff_exchange.pending_depth) == pointer(O.interior(state.runoff_pending_depth))
    O.TimeSteppers.time_step!(router,900)
    @test O.prognostic_state(router) == expected
    println("PUBLIC_NATIVE_RUNOFF_RESULT device=$(USE_GPU ? "gpu" : "cpu") native_depth=$(extrema(native)) endpoint_depth=$(extrema(endpoint)) exported_kg=$exported full_coupled=false")
end
println("PUBLIC_NATIVE_RUNOFF_PASS device=$(USE_GPU ? "gpu" : "cpu") full_coupled=false")
end
