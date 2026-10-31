using Dates, NCDatasets
include("test_sample.jl")
const NC = NCDatasets
const W = WriterIO.WC

@testset "Extended fields NetCDF roundtrip" begin
    weights=[2. 0.; 3. 5.]
    g=WriterIO.WriterGeometry.geometry([0.,180.],[-45.,45.],[-2.,-1.],
        [2. 7.; 3. 5.],Int8[1 0; 1 1],[0.,180.],[-45.,45.],[.25,.75])
    fields,scalars=additional_samples(API,fixture(),weights,[.25,.75]; ocean_latitude=LATITUDE)
    initial,initial_scalars=additional_samples(API,fixture(solved=false,iteration=0),weights,[.25,.75]; ocean_latitude=LATITUDE)
    annual=Dict("ocean_temperature"=>zeros(2,2,2),"ocean_salinity"=>fill(35.,2,2,2))
    cal=W.Calendar(DateTime(1993,1,1),3600)
    state=W.State(cal)
    for step in 1:24
        W.record!(state,step,fields,scalars,annual)
    end
    mktempdir() do directory
        launch=WriterIO.write_initial(NC,directory,cal,initial,initial_scalars,Dict();geometry=g)
        NC.NCDataset(launch) do ds
            @test all(isnan,ds["toa_net_downward"][:])
            @test ds["ice_volume_per_area"][:,:]==initial["ice_volume_per_area"]
        end
        events=W.finalize!(state)
        for event in events
            path=WriterIO.write_event(NC,directory,event,cal,Dict();geometry=g)
            event.kind==:month || continue
            NC.NCDataset(path) do ds
                for (name,values) in fields
                    @test Array(ds[name][:])==vec(values)
                    @test ds[name].attrib["units"]==LongRunFields.SPECS[name][2]
                end
                @test ds["daily_global_ice_volume"][:]==[5.]
                @test ds["daily_global_toa_net_downward"][:]==[10.]
                @test ds.attrib["sample_count"]==24
                @test ds.attrib["sampling_contract"]==WriterIO.CONTRACT
            end
        end
    end
end
println("LONG_RUN_FIELDS_IO_PASS")
