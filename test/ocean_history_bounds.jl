@testset "Ocean physical bounds include recovered historical excursions" begin
    check = ReadyESM._validate_ocean_extrema_history
    times = [0.0, 1.0, 2.0]
    Tmin, Tmax, Smin, Smax = fill(-2.0,3), fill(30.0,3), fill(20.0,3), fill(40.0,3)
    @test isnothing(check(times,Tmin,Tmax,Smin,Smax))
    for FT in (Float32, Float64)
        @test isnothing(check(FT.(times),fill(FT(-5),3),fill(FT(40),3),zeros(FT,3),fill(FT(50),3)))
        for (which, bad) in ((1,-5.1), (2,40.1), (3,-0.01), (4,50.1))
            fields = [FT.(x) for x in (Tmin,Tmax,Smin,Smax)]
            fields[which][2] = FT(bad) # Recovered endpoint must not hide this.
            @test_throws ErrorException check(times,fields...)
        end
        for bad in (NaN, Inf, -Inf), which in 1:4
            fields = [FT.(x) for x in (Tmin,Tmax,Smin,Smax)]
            fields[which][2] = FT(bad)
            @test_throws ErrorException check(times,fields...)
        end
    end
    # Reproduce the day363-to390 candidate's endpoint-only false pass.
    bad_s = [-0.1273276955, -0.2085178643, -0.2392190844, 0.4272964001]
    err = try
        check([363.,364.,365.,390.],fill(-2.,4),fill(36.,4),bad_s,fill(44.,4))
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("day 363.0",sprint(showerror,err))
    @test bad_s == [-0.1273276955, -0.2085178643, -0.2392190844, 0.4272964001]
    @test_throws ErrorException check(Float64[],Float64[],Float64[],Float64[],Float64[])
    @test_throws ErrorException check(times,Tmin[1:2],Tmax,Smin,Smax)
    @test_throws ErrorException check([0.,1.,1.],Tmin,Tmax,Smin,Smax)
    @test_throws ErrorException check([0.,2.,1.],Tmin,Tmax,Smin,Smax)
    @test_throws ErrorException check([0.,NaN,2.],Tmin,Tmax,Smin,Smax)
    @test_throws ErrorException check(times,fill(31.,3),Tmax,Smin,Smax)
    @test_throws ErrorException check(times,Tmin,Tmax,fill(41.,3),Smax)
    source = read(joinpath(ReadyESM.PROJECT_ROOT,"src","dynamic_ocean.jl"),String)
    @test occursin("_validate_ocean_extrema_history(\n        diagnostics.coupled_budget_time_days",source)
end
