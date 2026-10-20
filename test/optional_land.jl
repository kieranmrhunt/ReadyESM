using Dates

@testset "Wet atmosphere with optional land" begin
    R = ReadyESM
    S = R.SpeedyWeather
    # The callback also runs after atmospheric initialization.
    @test R._initialize_terrarium_temperature_from_atmosphere!(
        (; prognostic=(;)), nothing) === nothing
    @test R._initialize_terrarium_temperature_from_atmosphere!(
        (; prognostic=(; land=(;))), nothing) === nothing

    for F in (Float32, Float64), prescribed_land in (false, true)
        sg = S.SpectralGrid(; NF=F, trunc=7, nlayers=8)
        land = prescribed_land ?
            S.DryLandModel(sg; temperature=S.ConstantLandTemperature(sg; mask=false)) :
            nothing
        model = S.PrimitiveWetModel(sg;
            dynamics=false, parameterizations=(), land,
            ocean=S.PrescribedOcean(sg), sea_ice=S.PrescribedSeaIce(sg),
            output=nothing, feedback=nothing)
        simulation = S.initialize!(model)
        for lf in (1, 2)
            S.set!(simulation; lf, temperature=F(280), humidity=F(0.001),
                pressure=log(F(100000)), u=0, v=0)
        end
        S.transform!(simulation.variables, 1, model; initialize=true)
        @test haskey(simulation.variables.prognostic, :land) == prescribed_land
        @test R._initialize_terrarium_temperature_from_atmosphere!(
            simulation.variables, model) === nothing

        for internal_sync in ("false", "true")
            wrapped = deepcopy(simulation.variables)
            native = deepcopy(simulation.variables)
            withenv("READYESM_TRACE_SPEEDY_PHASES"=>"false",
                    "READYESM_SYNC_SPEEDY_INTERNAL_PHASES"=>internal_sync,
                    "READYESM_SYNC_SPEEDY_HANDOFF"=>"false") do
                for step in 1:4
                    dt = model.time_stepping.Δt
                    interval = step==1 ? dt/2 : step==2 ? dt : 2dt
                    lf1, lf2 = step<=2 ? 1 : 2, step==1 ? 1 : 2
                    S.timestep!(wrapped, interval, model, lf1, lf2)
                    invoke(S.timestep!,
                        Tuple{S.Variables,Real,S.PrimitiveEquation,Integer,Integer},
                        native, interval, model, lf1, lf2)
                    for name in (:vorticity, :divergence, :temperature, :humidity, :pressure)
                        @test isequal(getproperty(wrapped.prognostic,name).data,
                                      getproperty(native.prognostic,name).data)
                    end
                    for name in (:temperature, :humidity, :pressure, :u, :v,
                                 :temperature_prev, :humidity_prev, :pressure_prev)
                        @test isequal(getproperty(wrapped.grid,name).data,
                                      getproperty(native.grid,name).data)
                    end
                end
            end
        end
    end
end
println("OPTIONAL_LAND_TEST_PASS")
