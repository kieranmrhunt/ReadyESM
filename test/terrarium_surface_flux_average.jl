module TerrariumSurfaceFluxAverageTests
using ReadyESM, Test, Dates
const R = ReadyESM
const T, O = R.Terrarium, R.Oceananigans
const GPU_TEST = get(ENV,"READYESM_LAND_FLUX_GPU","false") == "true"
GPU_TEST && R.CUDA.allowscalar(false)
host(x) = Array(O.interior(x))
device(x) = GPU_TEST ? R.CUDA.CuArray(x) : x

function land(;FT=Float64,temp=20.,humidity=0.005,saturation=0.8)
    arch = GPU_TEST ? O.GPU() : O.CPU()
    grid = T.ColumnGrid(arch,FT,T.UniformSpacing(Δz=0.25,N=4))
    soil = T.SoilEnergyWaterCarbon(FT;hydrology=R._prognostic_soil_hydrology(FT))
    evap = R.ERA5PrescribedVegetationEvapotranspiration(FT)
    surface_hydrology = T.SurfaceHydrology(FT;evapotranspiration=evap,
        canopy_interception=T.NoCanopyInterception(FT))
    model = T.LandModel(;grid,soil,surface_hydrology,vegetation=nothing,snow=nothing)
    state = T.StateVariables(model)
    init = (;temperature=temp,saturation_water_ice=saturation,
        air_temperature=temp,specific_humidity=humidity,air_pressure=101325.,
        windspeed=2.,rainfall=0.,snowfall=0.,surface_shortwave_down=300.,
        surface_longwave_down=350.,vegetation_fraction=0.8,leaf_area_index=3.,
        root_fraction=reshape(FT[0.1,0.2,0.3,0.4],1,1,4))
    integrator = T.ModelIntegrator(state.clock,model,T.InputSources(FT),state,init)
    T.initialize!(integrator)
    return integrator
end

function reference!(integrator, seconds, step)
    schedule = R._terrarium_substep_schedule(Second(seconds),step)
    sums = Dict(:sensible_heat_flux=>0.,:latent_heat_flux=>0.)
    steps = fill(step,schedule.nfull)
    iszero(schedule.remainder) || push!(steps,schedule.remainder)
    for dt in steps
        T.timestep!(integrator,dt;finalize=false)
        for name in keys(sums)
            sums[name] += dt*Float64(only(host(getproperty(integrator.state,name))))
        end
    end
    T.compute_auxiliary!(integrator.state,integrator.model)
    return sums
end

@testset "Actual land substep fluxes and unchanged land trajectory" begin
    for FT in (Float32,Float64), (name,args) in
        (("warm",(;)),("dry",(;saturation=0.2)),
         ("condensing",(;humidity=0.03)),("frozen",(;temp=-5.,humidity=0.002)))
        for (seconds,dt) in ((60,60.),(75,60.),(900,60.),(900,6.))
            a,b = land(;FT,args...),land(;FT,args...)
            sums = reference!(a,seconds,dt)
            fluxes = R._terrarium_surface_flux_accumulator(b.state)
            # Reusing a scratch object must reset all previous interval totals.
            for buf in values(fluxes)
                fill!(buf,FT(1e6))
            end
            R._run_terrarium_land!(b,Second(seconds),dt;surface_fluxes=fluxes)
            @test isequal(map(host,a.state.prognostic),map(host,b.state.prognostic))
            @test a.state.clock.time == b.state.clock.time
            for field in (:sensible_heat_flux,:latent_heat_flux,:skin_temperature,:ground_heat_flux)
                @test isequal(host(getproperty(a.state,field)),host(getproperty(b.state,field)))
            end
            for name in keys(sums)
                applied = Float64(only(Array(getproperty(fluxes,name))))*seconds
                @test applied ≈ sums[name] rtol=512eps(FT) atol=512eps(FT)
            end
            if FT == Float64 && name=="warm" && seconds==900 && dt==60
                old_H = Float64(only(host(b.state.sensible_heat_flux)))*seconds
                old_LE = Float64(only(host(b.state.latent_heat_flux)))*seconds
                @test abs(old_H-sums[:sensible_heat_flux]) > 1e-3
                @test abs(old_LE-sums[:latent_heat_flux]) > 1e-3
                println("LAND_FLUX_AUDIT old_H_error_Wm2=",(old_H-sums[:sensible_heat_flux])/seconds,
                    " old_LE_error_Wm2=",(old_LE-sums[:latent_heat_flux])/seconds,
                    " new_H_error_Wm2=",only(Array(fluxes.sensible_heat_flux))-sums[:sensible_heat_flux]/seconds,
                    " new_LE_error_Wm2=",only(Array(fluxes.latent_heat_flux))-sums[:latent_heat_flux]/seconds)
            end
        end
        println("LAND_FLUX_CASE_PASS ",FT," ",name)
    end
end

@testset "Scratch reuse, duration rejection and array scatter" begin
    a,b = land(),land()
    fluxes = R._terrarium_surface_flux_accumulator(b.state)
    @test_throws ArgumentError R._run_terrarium_land!(b,Second(0),60.;surface_fluxes=fluxes)
    @test_throws ArgumentError R._run_terrarium_land!(b,Second(60),0.;surface_fluxes=fluxes)
    @test_throws ArgumentError R._average_terrarium_surface_fluxes!(fluxes,NaN)
    for seconds in (75,135)
        sums = reference!(a,seconds,60.)
        R._run_terrarium_land!(b,Second(seconds),60.;surface_fluxes=fluxes)
        @test isequal(map(host,a.state.prognostic),map(host,b.state.prognostic))
        for name in keys(sums)
            @test Float64(only(Array(getproperty(fluxes,name))))*seconds ≈ sums[name] rtol=1e-12 atol=1e-10
        end
    end
    destination = device(fill(-1.,6))
    source = device(reshape([2.,3.],2,1,1))
    indices = device([2,5])
    backend = R.KernelAbstractions.get_backend(destination)
    R._scatter_terrarium_surface_values!(backend,destination,source,indices;multiplier=2.,offset=1.)
    R.KernelAbstractions.synchronize(backend)
    @test Array(destination) == [-1.,5.,-1.,-1.,7.,-1.]
end
println("LAND_FLUX_AVERAGE_PASS gpu=",GPU_TEST)
end
