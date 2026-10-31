using Test, Dates, Serialization
const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(PROJECT_ROOT, "scripts", "climate_writer", "profiles", "profiles.jl"))
const P = ClimateProfiles
include(joinpath(PROJECT_ROOT, "scripts", "climate_writer", "monthly_writer", "writer_core.jl"))
const W = WriterCore

include("fixture.jl")

@testset "Pressure sampling, finite coverage and unchanged inputs" begin
    for FT in (Float64,Float32)
        args = fixture(;FT)
        before = deepcopy(args)
        fields = P.samples(args...)
        @test isequal(before,args)
        @test Set(keys(fields)) == Set(keys(P.SPECS))
        @test length(fields) == 12
        @test all(v -> all(isfinite,v),values(fields))
        @test fields["column_water"] ≈ Float64.(args[4]) .* 0.004 ./ 10 rtol=1e-7
        @test fields["surface_pressure"] == Float64.(args[4])
        for (name, kind, level, _) in P.TARGETS, i in 1:3
            valid = Float64(args[5][1])*Float64(args[4][i]) <= level <= Float64(args[5][end])*Float64(args[4][i])
            @test fields[name*"_valid_fraction"][i] == Float64(valid)
            if valid
                expected = kind == :temperature ? 250 + 3*log(level) : -80 + 7*log(level)
                @test fields[name*"_contribution"][i] ≈ expected atol=(FT==Float64 ? 1e-12 : 2e-5)
                @test P.complete_mean(fields,name)[i] == fields[name*"_contribution"][i]
            else
                @test fields[name*"_contribution"][i] == 0
                @test isnan(P.complete_mean(fields,name)[i])
            end
        end
        zero_args=deepcopy(args)
        fill!(zero_args[2],0); fill!(zero_args[3],0)
        zerosampled=P.samples(zero_args...)
        @test all(iszero,zerosampled["column_water"])
        @test zerosampled["zonal_wind_200hpa_contribution"] == zeros(3)
        @test zerosampled["zonal_wind_200hpa_valid_fraction"] == fields["zonal_wind_200hpa_valid_fraction"]
    end
    args=fixture()
    @test_throws DimensionMismatch P.samples(args[1],args[2][:,1:3],args[3:end]...)
    for index in (1,2,3,4,5,6)
        bad=deepcopy(args); bad[index][1]=NaN
        @test_throws ErrorException P.samples(bad...)
    end
    @test_throws ErrorException P.samples(args[1:6]...,0.)
    bad=deepcopy(args); bad[5][2]=bad[5][1]
    @test_throws ErrorException P.samples(bad...)
    bad=deepcopy(args); bad[6].*=.5
    @test_throws ErrorException P.samples(bad...)
    bad=deepcopy(args); bad[4][1]=0
    @test_throws ErrorException P.samples(bad...)
    reg=Dict{String,Tuple{Symbol,String,String}}()
    P.register_metadata!(reg)
    @test reg == P.SPECS
    @test_throws ErrorException P.register_metadata!(reg)
    @test_throws ArgumentError P.complete_mean(P.samples(args...),"unregistered")
    fields=P.samples(args...)
    fields["temperature_850hpa_valid_fraction"][1]=1.01
    @test_throws ErrorException P.complete_mean(fields,"temperature_850hpa")
end

@testset "Calendar means and split restart retain missing-sample coverage" begin
    cal=W.Calendar(DateTime(1992,12,31),21600)
    full=W.State(cal)
    split=W.State(cal)
    full_events=[]; split_events=[]
    for step in 1:8
        ps = step==2 ? [70000.,80000.,20000.] : [100000.,80000.,20000.]
        fields=P.samples(fixture(;pressure=ps)...)
        daily=Dict("dummy"=>[Float64(step)])
        annual=Dict("dummy"=>[Float64(step)])
        append!(full_events,W.record!(full,step,fields,daily,annual))
        append!(split_events,W.record!(split,step,fields,daily,annual))
        if step==3 || step==4
            buffer=IOBuffer(); serialize(buffer,W.snapshot(split)); seekstart(buffer)
            split=W.restore_state(deserialize(buffer),cal,step)
            @test isequal(W.snapshot(split),W.snapshot(full))
        end
    end
    @test isequal(full_events,split_events)
    @test isequal(W.snapshot(full),W.snapshot(split))
    december=only(e for e in full_events if e.kind==:month)
    @test december.key == (1992,12)
    @test december.count == 4 && !december.complete
    @test december.means["temperature_850hpa_valid_fraction"] == [.75,0,0]
    @test all(isnan,P.complete_mean(december.means,"temperature_850hpa"))
    @test P.complete_mean(december.means,"temperature_500hpa")[1] ≈ 250+3*log(50000.)
    @test split.month.count == 4
    @test split.month.sums["temperature_850hpa_valid_fraction"] == [4.,0,0]
end
println("CLIMATE_PROFILE_UNIT_PASS")
