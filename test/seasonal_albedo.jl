module SeasonalAlbedoTests
using ReadyESM, Test, Dates
const R=ReadyESM
const O=R.Oceananigans
fields(x)=(;(name=>getfield(x,name) for name in fieldnames(typeof(x)))...)
with(x;kwargs...)=R.ExperimentConfig(;merge(fields(x),(;kwargs...))...)

@testset "Spectral albedo configuration" begin
    released=R.load_config(joinpath(R.PROJECT_ROOT,"config/production.yml");check_input_files=false)
    @test released.terrarium_surface_albedo==:era5_spectral
    base=with(released;terrarium_surface_albedo=:constant,
        terrarium_albedo_cycle_path="",terrarium_albedo_cycle_sha256="")
    @test base.terrarium_surface_albedo==:constant
    @test isempty(R._seasonal_albedo_provenance(base))
    @test R.validate(base;check_input_files=false)===base
    good=with(base;terrarium_surface_albedo=:era5_spectral,
        terrarium_albedo_cycle_path="unused.nc",terrarium_albedo_cycle_sha256=repeat("a",64))
    @test R.validate(good;check_input_files=false)===good
    @test_throws ArgumentError R.validate(good)
    for change in ((;terrarium_surface_albedo=:unknown),
                   (;terrarium_albedo_cycle_path=""),
                   (;terrarium_albedo_cycle_sha256="bad"),
                   (;terrarium_albedo_crossing_band=:unknown),
                   (;land_model=:bucket))
        @test_throws ArgumentError R.validate(with(good;change...);check_input_files=false)
    end
    for rule in (:equal_halves,:visible,:near_infrared)
        @test R.validate(with(good;terrarium_albedo_crossing_band=rule);check_input_files=false) isa R.ExperimentConfig
    end
    for change in ((;terrarium_albedo_cycle_path="unused.nc"),
                   (;terrarium_albedo_cycle_sha256=repeat("a",64)),
                   (;terrarium_albedo_crossing_band=:visible))
        @test_throws ArgumentError R.validate(with(base;change...);check_input_files=false)
    end
    @test R._prepare_seasonal_albedo!(nothing,nothing,(;callbacks=Dict()))===nothing
    @test R._finish_seasonal_albedo_solve!(nothing,nothing,(;callbacks=Dict()))===nothing
    off=(;model=(;callbacks=Dict()))
    @test R._seasonal_albedo_restart_state(off)==(;)
    @test R._validate_seasonal_albedo_restart(off,(;))===nothing
    @test_throws ArgumentError R._validate_seasonal_albedo_restart(off,(;seasonal_albedo=(;)))
end

@testset "Held native input and default boundary" begin
    released=R.load_config(joinpath(R.PROJECT_ROOT,"config/production.yml");check_input_files=false)
    base=with(released;terrarium_surface_albedo=:constant,
        terrarium_albedo_cycle_path="",terrarium_albedo_cycle_sha256="")
    for FT in (Float32,Float64), RF in (Float32,Float64)
        old=R._coupled_land_surface_energy_balance(FT)
        off=R._coupled_land_surface_energy_balance(FT,base,RF)
        @test typeof(old)==typeof(off)
        @test old.albedo==off.albedo
        held=R.HeldSpectralLandAlbedo(FT,RF)
        @test held.initial_albedo==FT(RF(R.COUPLED_LAND_SURFACE_ALBEDO))
        @test held.emissivity==FT(RF(R.COUPLED_SURFACE_EMISSIVITY))
        grid=O.RectilinearGrid(O.CPU(),FT;size=(2,1,1),extent=(2,1,1),halo=(1,1,1))
        input=O.CenterField(grid)
        R._copy_held_albedo!(O.interior(input),RF[.2,.3])
        before=copy(parent(input))
        @test R.Terrarium.albedo(1,1,grid,(;albedo=input),held)==FT(RF(.2))
        @test R.Terrarium.compute_auxiliary!(nothing,grid,held)===nothing
        @test parent(input)==before
    end
end

@testset "Spectral boundary mapping and closure" begin
    A=R.SeasonalAlbedoAdapter
    bounds=[20000. 10000. 1000.; 25000. 18000. 10000.]
    @test A.visible_weights(bounds;crossing=:equal_halves)==[1.0,.5,0.0]
    @test A.visible_weights(bounds;crossing=:visible)==[1.0,1.0,0.0]
    @test A.visible_weights(bounds;crossing=:near_infrared)==[1.0,0.0,0.0]
    direct=zeros(3,2);diffuse=zeros(3,2)
    A.pack_albedos!(direct,diffuse,[1.0,.5,0.0],[.25,0.0],[.1,.2],
        [.6,NaN],[.8,NaN],[.4,NaN],[.2,NaN])
    @test direct[:,1]≈[.225,.25,.275]
    @test diffuse[:,1]≈[.175,.15,.125]
    @test direct[:,2]==fill(.2,3) && diffuse[:,2]==fill(.2,3)
    @test A.equivalent_land_albedo(Float32,.5,.1,100.,35.;previous=.4).albedo≈.6f0
    weak=A.equivalent_land_albedo(Float32,0.,.1,100.,10.;previous=.4)
    @test weak.albedo==.4f0 && weak.weak_signal
    @test_throws ArgumentError A.equivalent_land_albedo(Float32,.5,.1,100.,80.;previous=.4)
end

@testset "Seasonal cycle and land remapping" begin
    C=R.SeasonalVegetation
    @test C.bracket(DateTime(1993,1,17)).lower==1
    @test C.bracket(DateTime(1993,1,17)).upper==2
    @test C.bracket(DateTime(1993,12,31)).upper==1
    fields=Dict{Symbol,Array{Float64,3}}()
    for name in R.SeasonalLandRemapping.FIELD_NAMES
        fields[name]=fill(name==:lsm ? 1.0 : name in (:cvh,:cvl) ? .0 :
            name in (:lai_hv,:lai_lv) ? 1.0 : .3,3,4,12)
    end
    fields[:cvh].=0.6; fields[:cvl].=0.2
    fields[:lai_hv].=3.; fields[:lai_lv].=1.
    source=R.SeasonalLandRemapping.source([-90.,0.,90.],[0.,90.,180.,270.],fields;identity="synthetic")
    mapped=R.SeasonalLandRemapping.remap(source,[45.],[ -45.])
    @test mapped.cover==fill(.8,1,12)
    @test mapped.leaf_index≈fill(2.5,1,12)
    @test mapped.albedo[:aluvp]==fill(.3,1,12)
    @test !any(mapped.fallback) && mapped.source_identity=="synthetic"

    # Four land-conditioned weights can sum to one Float64 ULP above one.
    # The January production diagnostics require the mapped cover to remain
    # within the physical range without weakening their strict validator.
    fields[:cvh].=1.; fields[:cvl].=0.
    full=R.SeasonalLandRemapping.source([-90.,0.,90.],[0.,90.,180.,270.],fields)
    longitude,latitude=1.1050422658177905,-79.88410250277877
    @test sum(R.SeasonalLandRemapping.weights(full,longitude,latitude,1).weights)>1
    bounded=R.SeasonalLandRemapping.remap(full,[longitude],[latitude])
    @test bounded.cover==ones(1,12)
    @test all(isapprox.(bounded.leaf_index,fill(3.,1,12)))
    @test R.SeasonalLandRemapping.bounded_cover(1+eps(Float64))==1
    @test_throws ArgumentError R.SeasonalLandRemapping.bounded_cover(1+1e-8)
    @test_throws ArgumentError R.SeasonalLandRemapping.bounded_cover(-1e-8)
    @test_throws ArgumentError R.SeasonalLandRemapping.bounded_cover(NaN)
end

@testset "Enabled restart identity" begin
    identity=(;schema=1,source_sha256="a",mapping_sha256="b",crossing_band=:equal_halves)
    controller=R.SeasonalSpectralAlbedoController(nothing,[1],[.5],identity,Ref(false),
        zeros(Float32,2,1),zeros(Float32,2,1))
    simulation=(;model=(;callbacks=Dict(:seasonal_albedo=>controller)))
    record=(;identity,has_solved=false,direct=zeros(Float32,2,1),diffuse=zeros(Float32,2,1))
    @test R._validate_seasonal_albedo_restart(simulation,(;seasonal_albedo=record))===nothing
    @test_throws ArgumentError R._validate_seasonal_albedo_restart(simulation,(;))
    changed=(;identity=(;identity...,crossing_band=:visible),has_solved=false,
        direct=zeros(Float32,2,1),diffuse=zeros(Float32,2,1))
    @test_throws ArgumentError R._validate_seasonal_albedo_restart(simulation,(;seasonal_albedo=changed))
end
end
println("SEASONAL_ALBEDO_UNIT_PASS")
