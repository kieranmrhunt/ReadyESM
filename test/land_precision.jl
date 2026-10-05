module LandPrecisionTests
using ReadyESM, Test
const R = ReadyESM

function small_callback(::Type{FT}) where FT
    profiles = NamedTuple{R._ATMOSPHERE_PROCESS_PROFILE_ACCUMULATOR_FIELDS}(
        map(_ -> zeros(Float32, 1, 2), R._ATMOSPHERE_PROCESS_PROFILE_ACCUMULATOR_FIELDS))
    return R.GlobalAtmosphereDiagnosticsCallback{
        Float32,Vector{Float32},Vector{FT},Matrix{Float32},Nothing,FT,Vector{FT}}(;
        timestep_seconds=900f0, surface_sigma_thickness=0.1f0,
        point_weights=Float32[1],land_point_weights=Float32[1],land_column_weights=FT[1],
        column_cumulative_surface_water_vapor=zeros(Float32,1),
        column_cumulative_precipitation=zeros(Float32,1),
        column_cumulative_land_precipitation=zeros(FT,1),
        column_cumulative_land_evapotranspiration=zeros(FT,1),
        column_cumulative_land_surface_runoff=zeros(FT,1),
        soil_layer_thickness=FT[0.05,0.10],soil_porosity=FT(0.49),
        profiles...,radiation=nothing)
end

@testset "Explicit land precision and restart histories" begin
    @test R.ExperimentConfig().terrarium_precision == :float32
    @test R.validate(R.ExperimentConfig(terrarium_precision=:float64);check_input_files=false) isa R.ExperimentConfig
    @test_throws ArgumentError R.validate(R.ExperimentConfig(terrarium_precision=:float16);check_input_files=false)
    mktempdir() do tmp
        file=joinpath(tmp,"config.yml")
        write(file,"terrarium_precision: float64\n")
        @test R.load_config(file;check_input_files=false).terrarium_precision == :float64
    end
    for FT in (Float32,Float64)
        cb=small_callback(FT)
        @test eltype(cb.time_days) == Float32
        @test eltype(cb.point_weights) == Float32
        @test eltype(cb.column_cumulative_precipitation) == Float32
        for name in (:land_column_weights,:soil_layer_thickness,
                     :column_cumulative_land_precipitation,
                     :column_cumulative_land_evapotranspiration,
                     :column_cumulative_land_surface_runoff,
                     :land_total_water_storage,:land_water_budget_residual,
                     :cumulative_land_precipitation,:cumulative_land_evapotranspiration,
                     :cumulative_land_surface_runoff,:cumulative_balanced_launch_land_water_source)
            @test eltype(getproperty(cb,name)) == FT
        end
        @test typeof(cb.soil_porosity) == FT
        append!(cb.land_total_water_storage,FT[8320,8320+1e-5])
        cb.column_cumulative_land_precipitation[1]=FT(1e-5)
        append!(cb.cumulative_land_precipitation,FT[0,1e-5])
        cb.balanced_launch_land_water_source=FT(1e-9)
        restored=small_callback(FT)
        saved=R._atmosphere_diagnostics_callback_state(cb)
        R._restore_atmosphere_diagnostics_callback!(restored,saved)
        @test R._atmosphere_diagnostics_callback_state(restored) == saved
        @test typeof(restored.balanced_launch_land_water_source) == FT
        if FT==Float64
            @test cb.land_total_water_storage[2] != cb.land_total_water_storage[1]
            @test abs(diff(cb.land_total_water_storage)[1]-1e-5) < 1e-11
        end
        storage=R._column_soil_water_storage(fill(FT(0.5),1,1,2),zeros(FT,1,1,1),
                                           cb.soil_layer_thickness,cb.soil_porosity)
        @test eltype(storage) == FT
        @test only(storage) ≈ FT(1000)*FT(0.49)*FT(0.5)*sum(cb.soil_layer_thickness)
    end
    # No implicit migration of old Float32 land checkpoints is added here.
    @test_throws ErrorException R._copy_terrarium_restart_entry!(
        R.Terrarium.EmptyCache{Float64}(),R.Terrarium.EmptyCache{Float32}(),"precision_mismatch")
end
println("LAND_PRECISION_UNIT_PASS")
end
