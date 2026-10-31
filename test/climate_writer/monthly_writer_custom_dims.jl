using Test, Dates, NCDatasets
const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(PROJECT_ROOT, "scripts", "climate_writer", "monthly_writer", "writer_io.jl"))
using .WriterIO
const G = WriterIO.WriterGeometry
const W = WriterIO.WC
const NC = NCDatasets

function fixture()
    return G.geometry(Float32[0,90,180],Float32[-45,45],Float32[-50,-5],fill(2f6,3,2),ones(Int8,3,2),
        Float32[0,180,0,180],Float32[-45,-45,45,45],fill(.25f0,4))
end

@testset "Writer fields with their own dimensions" begin
    g=fixture()
    G.SPECS["test_psi"]=(:custom,"kg s-1","test streamfunction")
    G.SPECS["test_transport"]=(:custom,"W","test transport")
    G.SPECS["test_native_v"]=(:custom,"m s-1","test native velocity")
    G.CUSTOM_DIMS["test_psi"]=[("test_basin",[1.,2.],"1"),("test_interface",[1.,2.,3.],"1"),
                               ("test_lat",[-10.,0.,10.,20.],"degrees_north")]
    G.CUSTOM_DIMS["test_transport"]=[("test_basin",[1.,2.],"1"),("test_lat",[-10.,0.,10.,20.],"degrees_north")]
    G.CUSTOM_DIMS["test_native_v"]=[("test_i",[1.,2.,3.],"1"),("test_j",[1.,2.,3.],"1"),("test_k",[1.,2.],"1")]
    G.STORAGE["test_native_v"]=Float32
    psi=reshape(Float64.(1:24),2,3,4) .* 1.0e9
    transport=reshape(Float64.(1:8),2,4) .* 1.0e14
    native_v=fill(0.1,3,3,2)                                   # one more j row than the tracer grid
    monthly=Dict("sst"=>reshape(Float32.(1:6),3,2),"test_psi"=>psi,"test_transport"=>transport)
    annual=Dict("ocean_temperature"=>fill(8f0,3,2,2),"test_native_v"=>native_v)
    cal=W.Calendar(DateTime(1993,1,1),3600)
    state=W.State(cal)
    for step in 1:24; W.record!(state,step,monthly,Dict("global_sst"=>[3.5]),annual); end
    events=W.finalize!(state)
    mktempdir() do dir
        for event in events
            path=write_event(NC,dir,event,cal,Dict();geometry=g)
            NC.NCDataset(path) do ds
                if event.kind==:month
                    @test NC.dimnames(ds["test_psi"])==("test_basin","test_interface","test_lat")
                    @test ds["test_psi"][:,:,:]==psi
                    @test ds["test_transport"][:,:]==transport
                    @test ds["test_lat"][:]==[-10.,0.,10.,20.]
                    @test ds["test_lat"].attrib["units"]=="degrees_north"
                    @test ds["test_psi"].attrib["units"]=="kg s-1"
                    @test eltype(ds["test_psi"].var)==Float64
                    @test eltype(ds["sst"].var)==Float64            # existing fields unchanged
                    @test ds["sst"][:,:]==Float64.(monthly["sst"])
                else
                    @test NC.dimnames(ds["test_native_v"])==("test_i","test_j","test_k")
                    @test eltype(ds["test_native_v"].var)==Float32  # per-field storage type
                    @test ds["test_native_v"][:,:,:]≈Float32.(native_v)
                    @test eltype(ds["ocean_temperature"].var)==Float64
                end
            end
        end
        # A dimension already defined with other coordinate values is refused, and so is a wrong shape.
        G.CUSTOM_DIMS["test_transport"]=[("test_basin",[1.,2.],"1"),("test_lat",[-10.,0.,10.,30.],"degrees_north")]
        bad=joinpath(dir,"bad"); mkpath(bad)
        @test_throws ErrorException write_event(NC,bad,first(e for e in events if e.kind==:month),cal,Dict();geometry=g)
        G.CUSTOM_DIMS["test_transport"]=[("test_basin",[1.,2.],"1"),("test_lat",[-10.,0.,10.],"degrees_north")]
        bad2=joinpath(dir,"bad2"); mkpath(bad2)
        @test_throws DimensionMismatch write_event(NC,bad2,first(e for e in events if e.kind==:month),cal,Dict();geometry=g)
        delete!(G.CUSTOM_DIMS,"test_transport")
        bad3=joinpath(dir,"bad3"); mkpath(bad3)
        @test_throws ErrorException write_event(NC,bad3,first(e for e in events if e.kind==:month),cal,Dict();geometry=g)
    end
end
