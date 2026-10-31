using Test, Dates, NCDatasets, Serialization
const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(PROJECT_ROOT, "scripts", "climate_writer", "monthly_writer", "writer_io.jl"))
include(joinpath(PROJECT_ROOT, "scripts", "climate_writer", "profiles", "sample.jl"))
include("fixture.jl")
const W=WriterIO.WC
const P=ClimateProfiles
const NC=NCDatasets

@testset "Adapter copies state without mutation" begin
    args=fixture()
    t,q,u,ps,sigma,thick,g=args
    atmosphere=(;variables=(;grid=(;temperature=(;data=t),humidity=(;data=q),u=(;data=u),pressure=(;data=log.(ps)))),
        model=(;geometry=(;σ_levels_full=sigma,σ_levels_thick=thick),planet=(;gravity=g)))
    before=deepcopy(atmosphere)
    actual=climate_profile_samples(nothing,atmosphere)
    @test isequal(before,atmosphere)
    expected=P.samples(t,q,u,exp.(log.(ps)),sigma,thick,g)
    @test actual==expected
    actual["surface_pressure"][1]=0
    @test isequal(before,atmosphere)
end

@testset "Complete monthly NetCDF contributions and coverage" begin
    output=ENV["CLIMATE_PROFILE_TEST_WORK"]
    ispath(output) && error("fresh test output required")
    mkpath(output)
    geometry=WriterIO.WriterGeometry.geometry([0.],[0.],[-10.],ones(1,1),ones(Int8,1,1),
        [0.,120.,240.],[-45.,0.,45.],[.25,.5,.25])
    cal=W.Calendar(DateTime(1993,1,1),21600)
    full,split=W.State(cal),W.State(cal)
    full_events=[];split_events=[]
    initial=P.samples(fixture()...)
    WriterIO.write_initial(NC,output,cal,initial,Dict{String,Any}(),Dict();geometry)
    for step in 1:124
        ps=step==1 ? [70000.,80000.,20000.] : [100000.,80000.,20000.]
        fields=P.samples(fixture(;pressure=ps)...)
        daily=Dict{String,Any}()
        append!(full_events,W.record!(full,step,fields,daily,fields))
        append!(split_events,W.record!(split,step,fields,daily,fields))
        if step==50
            buffer=IOBuffer();serialize(buffer,W.snapshot(split));seekstart(buffer)
            split=W.restore_state(deserialize(buffer),cal,step)
        end
    end
    @test isequal(full_events,split_events)
    @test isequal(W.snapshot(full),W.snapshot(split))
    event=only(full_events)
    @test event.kind==:month && event.complete && event.count==124
    @test event.start==DateTime(1993,1,1) && event.stop==DateTime(1993,2,1)
    metadata=Dict("profile_schema"=>"climate_pressure_contributions_v1",
        "profile_interpolation"=>"linear in log pressure, no extrapolation",
        "profile_coverage"=>"contribution and coverage averaged over every coupled step; score only valid_fraction == 1")
    file=WriterIO.write_event(NC,output,event,cal,metadata;geometry)
    NC.NCDataset(file) do ds
        @test ds.attrib["complete"]==1 && ds.attrib["sample_count"]==124
        @test ds.attrib["sampling_contract"]==WriterIO.CONTRACT
        @test ds.attrib["profile_schema"]=="climate_pressure_contributions_v1"
        for (name,values) in event.means
            @test Array(ds[name][:])==values
            @test ds[name].attrib["units"]==P.SPECS[name][2]
            @test ds[name].attrib["coordinates"]=="atmosphere_longitude atmosphere_latitude"
        end
        fields=Dict(name=>Array(ds[name][:]) for name in keys(P.SPECS))
        @test fields["temperature_850hpa_valid_fraction"]==[123/124,0,0]
        @test all(isnan,P.complete_mean(fields,"temperature_850hpa"))
        @test P.complete_mean(fields,"temperature_500hpa")[1] ≈ 250+3*log(50000.)
        @test fields["column_water"][1] ≈ (123*40+28)/124
        @test ds["day_count"][:]==fill(4,31)
    end
end
println("CLIMATE_PROFILE_IO_PASS")
