using Test, Dates, NCDatasets
const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(PROJECT_ROOT, "scripts", "climate_writer", "monthly_writer", "writer_io.jl"))
using .WriterIO
const G = WriterIO.WriterGeometry
const W = WriterIO.WC
const NC = NCDatasets

function fixture(;geographic=true)
    x=Float32[0,90,180]
    y=Float32[-45,45]
    z=Float32[-50,-5]
    area=fill(2f6,3,2)
    wet=ones(Int8,3,2); wet[1,1]=0
    lon=Float32[0,180,0,180]
    lat=Float32[-45,-45,45,45]
    weights=fill(.25f0,4)
    return G.geometry(x,y,z,area,wet,lon,lat,weights;geographic)
end

@testset "Writer native-grid metadata" begin
    g=fixture()
    @test g.x[:,1] == [0,90,180]
    @test g.y[1,:] == [-45,45]
    @test sum(g.weights)==1
    @test g.active[1,1]==0
    for (field,value) in ((:area,fill(-1.,3,2)),(:active,fill(2,3,2)),
                         (:z,[0.,-5.]),(:latitude,[91.,0.,0.,0.]),(:weights,ones(4)))
        args=merge(g,NamedTuple{(field,)}((value,)))
        @test_throws ErrorException G.geometry(args.x,args.y,args.z,args.area,args.active,
            args.longitude,args.latitude,args.weights)
    end
    @test_throws DimensionMismatch G.geometry(zeros(2),zeros(2),[-1.],ones(3,2),ones(3,2),[0.],[0.],[1.])

    monthly=Dict("sst"=>reshape(Float32.(1:6),3,2),"sss"=>fill(35f0,3,2),
        "ice_concentration"=>fill(.5f0,3,2),"air_temperature_lowest"=>fill(280f0,4),
        "outgoing_longwave"=>fill(240f0,4),"rain_rate"=>fill(1f-8,4))
    annual=Dict("ocean_temperature"=>fill(8f0,3,2,2),"ocean_salinity"=>fill(35f0,3,2,2))
    daily=Dict("global_sst"=>[3.5],"global_outgoing_longwave"=>[240.])
    cal=W.Calendar(DateTime(1993,1,1),3600)
    state=W.State(cal)
    for step in 1:24; W.record!(state,step,monthly,daily,annual); end
    events=W.finalize!(state)
    mktempdir() do dir
        old=joinpath(dir,"old"); new=joinpath(dir,"new")
        mkpath(old); mkpath(new)
        a=write_initial(NC,old,cal,monthly,daily,Dict())
        b=write_initial(NC,new,cal,monthly,daily,Dict();geometry=g)
        NC.NCDataset(a) do da
            NC.NCDataset(b) do db
                for name in keys(monthly); @test da[name][:] == db[name][:]; end
                @test db["rain_rate"].attrib["units"]=="m s-1"
                @test db["sst"].attrib["units"]=="degree_Celsius"
                @test db["air_temperature_lowest"].attrib["units"]=="K"
                @test db["sst"].attrib["sampling_statistic"]=="instantaneous launch state"
                @test db["ocean_latitude"][:,:]==g.y
                @test db["ocean_active_mask"][:,:]==g.active
                @test db["atmosphere_area_weight"][:]==g.weights
                @test db["scalar_global_sst"][]==3.5
                @test db["scalar_global_sst"].attrib["units"]=="degree_Celsius"
                @test NC.dimnames(db["sst"])==("ocean_i","ocean_j")
            end
        end
        for event in events
            a=write_event(NC,old,event,cal,Dict())
            b=write_event(NC,new,event,cal,Dict();geometry=g)
            NC.NCDataset(a) do da
                NC.NCDataset(b) do db
                    for name in keys(event.means); @test da[name][:] == db[name][:]; end
                    @test db.attrib["sample_count"]==24
                    @test db.attrib["complete"]==0
                    @test db["ocean_z"][:]==[-50.,-5.]
                    @test db["ocean_z"].attrib["positive"]=="up"
                    if event.kind==:month
                        @test db["daily_global_sst"][:]==[3.5]
                        @test db["daily_global_outgoing_longwave"].attrib["units"]=="W m-2"
                    else
                        @test NC.dimnames(db["ocean_temperature"])==("ocean_i","ocean_j","ocean_k")
                        @test db["ocean_temperature"].attrib["coordinates"]=="ocean_longitude ocean_latitude ocean_z"
                    end
                end
            end
        end
        cart=joinpath(dir,"cartesian"); mkpath(cart)
        path=write_initial(NC,cart,cal,monthly,daily,Dict();geometry=fixture(geographic=false))
        NC.NCDataset(path) do ds
            @test haskey(ds,"ocean_x") && !haskey(ds,"ocean_longitude")
            @test ds["ocean_x"].attrib["units"]=="m"
            @test ds["sst"].attrib["coordinates"]=="ocean_x ocean_y"
        end
    end
end
println("MONTHLY_WRITER_GEOMETRY_PASS")
