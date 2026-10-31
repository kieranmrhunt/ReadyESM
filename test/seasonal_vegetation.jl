module SeasonalVegetationTests
using ReadyESM, Test, Dates
const R=ReadyESM
const O=R.Oceananigans
fields(x)=(;(name=>getfield(x,name) for name in fieldnames(typeof(x)))...)
with(x;kwargs...)=R.ExperimentConfig(;merge(fields(x),(;kwargs...))...)

@testset "Seasonal vegetation configuration" begin
    released=R.load_config(joinpath(R.PROJECT_ROOT,"config/production.yml");check_input_files=false)
    @test released.terrarium_vegetation_cycle==:cover_lai
    base=with(released;terrarium_vegetation_cycle=:off,
        terrarium_vegetation_cycle_path="",terrarium_vegetation_cycle_sha256="")
    @test base.terrarium_vegetation_cycle==:off
    @test R.validate(base;check_input_files=false)===base
    good=with(base;terrarium_vegetation_cycle=:cover_lai,
        terrarium_vegetation_cycle_path="unused.nc",terrarium_vegetation_cycle_sha256=repeat("a",64))
    @test R.validate(good;check_input_files=false)===good
    @test_throws ArgumentError R.validate(good)
    for change in ((;terrarium_vegetation_cycle=:unknown),
                   (;terrarium_vegetation_cycle_path=""),
                   (;terrarium_vegetation_cycle_sha256="bad"),
                   (;terrarium_evapotranspiration=:bare_ground))
        @test_throws ArgumentError R.validate(with(good;change...);check_input_files=false)
    end
    @test_throws ArgumentError R.validate(with(base;terrarium_vegetation_cycle_path="unused.nc");check_input_files=false)
end

@testset "One input update, unchanged root and stores" begin
    monthly=(;cover=fill(.5,2,12),area_leaf_index=fill(1.5,2,12),source_identity=repeat("a",64))
    monthly.area_leaf_index[:,7].=2.5
    for FT in (Float32,Float64)
        controller=R.SeasonalVegetationController(monthly,zeros(FT,2),zeros(FT,2),repeat("b",64))
        grid=O.RectilinearGrid(O.CPU(),FT;size=(2,1,1),extent=(2,1,1),halo=(1,1,1))
        cover=O.CenterField(grid);lai=O.CenterField(grid);root=O.CenterField(grid)
        water=O.CenterField(grid);energy=O.CenterField(grid)
        O.set!(root,.25);O.set!(water,.3);O.set!(energy,4)
        inputs=(;vegetation_fraction=cover,leaf_area_index=lai,root_fraction=root)
        before=copy.((parent(root),parent(water),parent(energy)))
        land=(;terrarium=(;inputs,prognostic=(;water,energy)))
        date=R.SeasonalVegetation.month_centre(DateTime(1993,1,1))
        vars=(;prognostic=(;land,clock=(;time=date)))
        model=(;callbacks=Dict(:seasonal_vegetation=>controller))
        R._update_seasonal_vegetation!(vars,model)
        @test all(==(.5),O.interior(cover))
        @test all(==(3),O.interior(lai))
        @test (parent(root),parent(water),parent(energy))==before
        summer=(;prognostic=(;land,clock=(;time=R.SeasonalVegetation.month_centre(DateTime(1993,7,1)))))
        R._update_seasonal_vegetation!(summer,model)
        @test all(==(5),O.interior(lai))
        @test (parent(root),parent(water),parent(energy))==before
        R._update_seasonal_vegetation!(vars,model)
        @test all(==(3),O.interior(lai))
        @test R._update_seasonal_vegetation!((;),(;callbacks=Dict()))===nothing
        sim=(;model)
        record=R._seasonal_vegetation_restart_state(sim)
        @test R._validate_seasonal_vegetation_restart(sim,record)===nothing
        @test_throws ArgumentError R._validate_seasonal_vegetation_restart(sim,(;))
        off=(;model=(;callbacks=Dict()))
        @test R._seasonal_vegetation_restart_state(off)==(;)
        @test R._validate_seasonal_vegetation_restart(off,(;))===nothing
        @test_throws ArgumentError R._validate_seasonal_vegetation_restart(off,record)
        wrong=(;seasonal_vegetation=merge(record.seasonal_vegetation,(;source_sha256=repeat("c",64))))
        @test_throws ArgumentError R._validate_seasonal_vegetation_restart(sim,wrong)
        wrong=(;seasonal_vegetation=merge(record.seasonal_vegetation,(;mapping_sha256=repeat("d",64))))
        @test_throws ArgumentError R._validate_seasonal_vegetation_restart(sim,wrong)
        @test (parent(root),parent(water),parent(energy))==before
    end
end
end
println("SEASONAL_VEGETATION_INPUT_UNIT_PASS")
