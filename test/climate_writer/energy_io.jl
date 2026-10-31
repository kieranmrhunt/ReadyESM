using Test, Dates, NCDatasets, Serialization
const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(PROJECT_ROOT, "scripts", "climate_writer", "monthly_writer", "writer_io.jl"))
include(joinpath(PROJECT_ROOT, "scripts", "climate_writer", "energy", "energy.jl"))
const E = ClimateEnergy
const W = WriterIO.WC
const NC = NCDatasets
E.register_metadata!(WriterIO.WriterGeometry.SPECS)

@testset "Monthly energy components and restart" begin
    output = ENV["CLIMATE_ENERGY_TEST_WORK"]
    ispath(output) && error("fresh output required")
    mkpath(output)
    geometry = WriterIO.WriterGeometry.geometry([0.,180.],[0.],[-10.],ones(2,1),ones(Int8,2,1),
        [0.,180.],[-45.,45.],[.5,.5])
    calendar = W.Calendar(DateTime(1993,1,1),21600)
    full, split = W.State(calendar), W.State(calendar)
    full_events, split_events = [], []
    expected = Dict{String,Any}()
    metadata = Dict("energy_schema"=>"native_latent_components_v1",
        "energy_ocean_timing"=>"latest post-step interface cache",
        "energy_land_timing"=>"retained land substep-average humidity flux times its vaporisation constant",
        "energy_initial_availability"=>"NaN, excluded from period means")
    initial = Dict{String,Any}(name=>fill(NaN,kind==:surface ? (2,1) : (2,))
        for (name,(kind,_,_)) in E.SPECS)
    WriterIO.write_initial(NC,output,calendar,initial,Dict{String,Any}(),metadata;geometry)
    for step in 1:124
        concentration = reshape([step/124,0.],2,1)
        fields = Dict{String,Any}(E.ocean_samples(fill(100.,2,1),fill(80.,2,1),concentration,ones(2,1)))
        merge!(fields,E.vapor_samples(fill(1e-5,2,1),fill(-2e-6,2,1),concentration,ones(2,1)))
        merge!(fields,E.land_samples([step*1e-6,NaN],[.25,0.],2.5e6))
        for (name,values) in fields
            get!(expected,name,zeros(size(values))) .+= values/124
        end
        daily = Dict{String,Any}()
        append!(full_events,W.record!(full,step,fields,daily,fields))
        append!(split_events,W.record!(split,step,fields,daily,fields))
        if step==50
            buffer=IOBuffer()
            serialize(buffer,W.snapshot(split))
            seekstart(buffer)
            split=W.restore_state(deserialize(buffer),calendar,step)
        end
    end
    @test isequal(full_events,split_events)
    @test isequal(W.snapshot(full),W.snapshot(split))
    event=only(full_events)
    @test event.kind==:month && event.complete && event.count==124
    @test event.start==DateTime(1993,1,1) && event.stop==DateTime(1993,2,1)
    path=WriterIO.write_event(NC,output,event,calendar,metadata;geometry)
    NC.NCDataset(path) do ds
        @test ds.attrib["complete"]==1 && ds.attrib["sample_count"]==124
        @test ds.attrib["energy_schema"]=="native_latent_components_v1"
        @test ds.attrib["sampling_contract"]==WriterIO.CONTRACT
        for (name,values) in expected
            @test Array(ds[name][ntuple(_ -> Colon(), ndims(ds[name]))...]) ≈ values
            @test ds[name].attrib["units"]==E.SPECS[name][2]
            @test ds[name].attrib["coordinates"] == (E.SPECS[name][1]==:surface ?
                "ocean_longitude ocean_latitude" : "atmosphere_longitude atmosphere_latitude")
        end
        @test ds["day_count"][:]==fill(4,31)
        @test ds["ocean_ice_latent_heat_upward"][:,:] ≈
            ds["open_water_latent_heat_contribution"][:,:]+ds["sea_ice_latent_heat_contribution"][:,:]
    end
end
println("CLIMATE_ENERGY_IO_TEST_PASS")
