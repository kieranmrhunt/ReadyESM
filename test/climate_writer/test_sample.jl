using Test
const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const WRITER_DIR = joinpath(PROJECT_ROOT, "scripts", "climate_writer", "monthly_writer")
include(joinpath(WRITER_DIR, "writer_io.jl"))
include(joinpath(WRITER_DIR, "writer_sampling.jl"))
include(joinpath(PROJECT_ROOT, "scripts", "climate_writer", "extended_fields", "sample.jl"))

module FakeOcean
interior(x) = x
end
module FakeRadiation
const reads = Ref(0)
for (fn, field) in ((:sw_flux_dn,:swd),(:sw_flux_up,:swu),(:lw_flux_dn,:lwd),
                    (:lw_flux_up,:lwu),(:clear_sw_flux_up,:cswu),(:clear_lw_flux_up,:clwu))
    @eval function $fn(solver)
        reads[] += 1
        return getproperty(solver, $(QuoteNode(field)))
    end
end
end

const LATITUDE = [-45. 45.; -45. 45.]
const API = (Oceananigans=FakeOcean, RRTMGP=FakeRadiation)
function fixture(; solved=true, iteration=1, clear=true)
    swd=Float32[100 120; 340 360]
    swu=Float32[10 12; 90 100]
    lwd=Float32[300 310; 0 0]
    lwu=Float32[390 400; 240 250]
    solver=(; swd,swu,lwd,lwu,cswu=swu.-5,clwu=lwu.+10,clear_flux_lw=clear ? true : nothing)
    radiation=(call_counter=solved ? 1 : 0,solver)
    p=(rain_rate=(data=Float32[1e-8,2e-8],),snow_rate=(data=Float32[3e-8,4e-8],),
       sensible_heat_flux=(data=Float32[-1,2],),surface_humidity_flux=(data=Float32[-2e-5,3e-5],))
    atmosphere=(model=(longwave_radiation=radiation,),variables=(parameterizations=p,))
    ice=(ice_thickness=reshape(Float32[2,4,6,8],2,2,1),
         ice_concentration=reshape(Float32[.5,.25,1,0],2,2,1))
    return (atmosphere=atmosphere,sea_ice=(model=ice,),clock=(iteration=iteration,))
end

@testset "Long-run sampler phase, mask and units" begin
    ocean_weights = [2. 0.; 3. 5.]
    atmosphere_weights = [.25,.75]
    earth=fixture()
    before=deepcopy(earth)
    FakeRadiation.reads[]=0
    monthly,daily=additional_samples(API,earth,ocean_weights,atmosphere_weights; ocean_latitude=LATITUDE)
    @test isequal(earth,before)
    @test FakeRadiation.reads[]==6
    @test daily["global_toa_net_downward"]==[10.]
    @test daily["global_ice_area"]==[1.75]
    @test daily["global_ice_volume"]==[5.]
    @test monthly["surface_sensible_heat_upward"]==[-1.,2.]
    @test monthly["total_precipitation_rate"]==Float64.(earth.atmosphere.variables.parameterizations.rain_rate.data).+
        Float64.(earth.atmosphere.variables.parameterizations.snow_rate.data)
    @test all(name -> haskey(WriterIO.WriterGeometry.SPECS,name), union(keys(monthly),keys(daily)))
    FakeRadiation.reads[]=0
    initial,initial_daily=additional_samples(API,fixture(solved=false,iteration=0),ocean_weights,atmosphere_weights; ocean_latitude=LATITUDE)
    @test FakeRadiation.reads[]==0
    @test all(isnan,initial["toa_net_downward"])
    @test isnan(only(initial_daily["global_toa_net_downward"]))
    @test initial["ice_volume_per_area"]==monthly["ice_volume_per_area"]
    @test_throws ErrorException additional_samples(API,fixture(solved=false),ocean_weights,atmosphere_weights; ocean_latitude=LATITUDE)
    @test_throws ErrorException additional_samples(API,fixture(clear=false),ocean_weights,atmosphere_weights; ocean_latitude=LATITUDE)
    @test daily["northern_ice_area"]==[0.]
    @test daily["southern_ice_area"]==[1.75]
    @test daily["northern_ice_volume"]==[0.]
    @test daily["southern_ice_volume"]==[5.]
    @test daily["northern_ice_area"] .+ daily["southern_ice_area"] == daily["global_ice_area"]
    @test daily["northern_ice_volume"] .+ daily["southern_ice_volume"] == daily["global_ice_volume"]
    both=fixture()
    both.sea_ice.model.ice_concentration[2,2,1]=.5f0
    _,both_daily=additional_samples(API,both,ocean_weights,atmosphere_weights; ocean_latitude=LATITUDE)
    @test both_daily["northern_ice_area"]==[2.5]
    @test both_daily["southern_ice_area"]==[1.75]
    @test both_daily["northern_ice_volume"]==[20.]
    @test both_daily["southern_ice_volume"]==[5.]
    _,equator=additional_samples(API,earth,ocean_weights,atmosphere_weights; ocean_latitude=zeros(2,2))
    @test equator["northern_ice_area"]==daily["global_ice_area"]
    @test equator["southern_ice_area"]==[0.]
end
println("LONG_RUN_SAMPLER_UNIT_PASS")
