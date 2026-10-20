module CarbonDiagnosticTests
using ReadyESM, Test, NCDatasets
const R = ReadyESM
const O = R.Oceananigans

function controller(; mode=:emissions_driven, rate=10.0, history=true)
    config = R.ExperimentConfig(; ocean_model=:oceananigans,
        forcing=R.ForcingConfig(radiation=:rrtmgp_clear_sky),
        carbon=R.CarbonConfig(; mode, emissions_pgc_per_year=rate))
    value = R.CarbonController(mode, 420.0 * R._CARBON_PGC_PER_PPM, 420.0,
        38000.0, 0.0, 0.0, 0.0, 0.0, NaN, -1, 900.0, R._carbon_provenance(config))
    history && R._initialize_carbon_diagnostic_history!(value)
    return value
end

function advance!(value, step; cadence=96, at_end=false)
    rate = 1e-6 * cos(step / 17)
    R._advance_carbon_ledger!(value, step, 900.0 * step, rate)
    R._record_carbon_diagnostics!(value; sample_every_n_steps=cadence, at_end)
    return value
end

function physical_state(value)
    state = O.prognostic_state(value)
    names = Tuple(name for name in keys(state) if name != :diagnostic_history)
    return NamedTuple{names}(Tuple(getproperty(state, name) for name in names))
end

@testset "Carbon history is passive, sampled and restartable" begin
    for (mode, rate) in ((:prescribed, 0.0), (:emissions_driven, 0.0),
                         (:emissions_driven, 10.0)), cadence in (1, 4, 96)
        live = controller(; mode, rate)
        bare = controller(; mode, rate, history=false)
        @test isequal(physical_state(live), physical_state(bare))
        @test isempty(R._carbon_diagnostic_history(live).rows)
        saved = nothing
        for step in 0:100
            advance!(live, step; cadence)
            advance!(bare, step; cadence)
            @test physical_state(live) == physical_state(bare)
            before = O.prognostic_state(live)
            R._record_carbon_diagnostics!(live; sample_every_n_steps=cadence)
            @test isequal(before, O.prognostic_state(live))
            step == 37 && (saved = O.prognostic_state(live))
        end
        history = R._carbon_diagnostic_history(live)
        @test [row.iteration for row in history.rows] == collect(0:cadence:100)
        before = O.prognostic_state(live)
        rows = R._carbon_output_rows(live)
        @test last(rows).iteration == 100
        @test isequal(before, O.prognostic_state(live))
        empty!(rows)
        @test !isempty(history.rows)
        @test last(saved.diagnostic_history.rows).iteration <= 37

        mktempdir() do directory
            file = joinpath(directory, "carbon.jld2")
            R.JLD2.jldsave(file; state=saved)
            for source in (saved, R.JLD2.load(file, "state"))
                restored = controller(; mode, rate)
                @test O.restore_prognostic_state!(restored, source) === restored
                @test isequal(O.prognostic_state(restored), saved)
                for step in 38:100
                    advance!(restored, step; cadence)
                end
                @test isequal(O.prognostic_state(restored), O.prognostic_state(live))
            end
        end
    end
end

@testset "Carbon history rejects incompatible checkpoints before mutation" begin
    c = controller()
    for step in 0:96
        advance!(c, step)
    end
    state = O.prognostic_state(c)
    for history in (
        merge(state.diagnostic_history, (; schema="unknown")),
        merge(state.diagnostic_history, (; interval_seconds=1.0)),
        merge(state.diagnostic_history, (; initial_total_pgc=NaN)),
        merge(state.diagnostic_history, (; rows=NamedTuple[])),
        merge(state.diagnostic_history, (; rows=reverse(state.diagnostic_history.rows))),
    )
        @test_throws ErrorException O.restore_prognostic_state!(c,
            merge(state, (; diagnostic_history=history)))
        @test isequal(O.prognostic_state(c), state)
    end
    @test_throws ErrorException O.restore_prognostic_state!(c, physical_state(c))
    @test isequal(O.prognostic_state(c), state)
    bare = controller(history=false)
    @test_throws ErrorException O.restore_prognostic_state!(bare, state)
    @test_throws ErrorException R._initialize_carbon_diagnostic_history!(c)
    @test_throws ErrorException R._record_carbon_diagnostics!(c; sample_every_n_steps=1)
    @test isequal(O.prognostic_state(c), state)
end

function output_fixture(mode)
    c = controller(; mode, rate=mode == :emissions_driven ? 10.0 : 0.0)
    rows = NamedTuple[]
    initial = nothing
    interval = nothing
    error = nothing
    if mode != :tracer_only
        for step in 0:100
            advance!(c, step)
        end
        rows = R._carbon_output_rows(c)
        initial = R._carbon_diagnostic_history(c).initial_total_pgc
        interval = 86400.0
        error = 1e-8
    end
    wet = trues(2, 2, 2)
    wet[1, 2, 1] = false
    dic = R._physical_carbon_output(fill(-1f0, 2, 2, 2), wet, 2000)
    ta = R._physical_carbon_output(fill(2f0, 2, 2, 2), wet, 2300)
    return (; carbon_diagnostics=(
        schema=R._CARBON_HISTORY_SCHEMA, mode=String(mode),
        endpoint_time_seconds=90000.0, interval_seconds=interval,
        initial_total_pgc=initial, ocean_tracer_inventory_pgc=38000.0,
        ocean_inventory_relative_error=error, reference_density_kg_m3=1026.0,
        dic_storage_offset_umol_kg=2000.0, alkalinity_storage_offset_umol_kg=2300.0,
        rows, wet, dissolved_inorganic_carbon_umol_kg=dic, alkalinity_umol_kg=ta))
end

@testset "Carbon maps, output endpoints and fixed NetCDF variables" begin
    @test R._carbon_diagnostic_fields(nothing, R.ExperimentConfig()) == NamedTuple()
    @test R._validate_carbon_diagnostics(NamedTuple()) === nothing
    source = fill(-1f0, 2, 2, 2)
    wet = trues(2, 2, 2)
    wet[1, 1, 1] = false
    physical = R._physical_carbon_output(source, wet, 2000)
    @test all(==(-1f0), source)
    @test all(==(1999.0), physical[wet])
    @test isnan(physical[1, 1, 1])
    @test_throws ErrorException R._physical_carbon_output(source, trues(1, 1, 1), 2000)

    variable_sets = Set{String}[]
    mktempdir() do directory
        for mode in (:prescribed, :emissions_driven, :tracer_only, :off)
            diagnostics = mode == :off ? NamedTuple() : output_fixture(mode)
            @test R._validate_carbon_diagnostics(diagnostics) === nothing
            file = joinpath(directory, "$mode.nc")
            NCDataset(file, "c") do ds
                for name in ("ocean_x", "ocean_y", "ocean_layer")
                    defDim(ds, name, 2)
                end
                R._write_carbon_diagnostics!(ds, diagnostics,
                    ("ocean_x", "ocean_y"), "1993-01-01 00:00:00")
            end
            NCDataset(file) do ds
                if mode == :off
                    @test isempty(keys(ds))
                    @test !haskey(ds.attrib, "carbon_mode")
                    return
                end
                @test ds.attrib["carbon_mode"] == String(mode)
                @test ds.attrib["carbon_dic_storage_offset_umol_kg"] == 2000
                @test ds.attrib["carbon_alkalinity_storage_offset_umol_kg"] == 2300
                @test isequal(ds["ocean_dissolved_inorganic_carbon"][:, :, :],
                              diagnostics.carbon_diagnostics.dissolved_inorganic_carbon_umol_kg)
                @test isequal(ds["ocean_total_alkalinity"][:, :, :],
                              diagnostics.carbon_diagnostics.alkalinity_umol_kg)
                if mode == :tracer_only
                    @test Set(keys(ds)) == Set(("ocean_dissolved_inorganic_carbon", "ocean_total_alkalinity"))
                else
                    push!(variable_sets, Set(keys(ds)))
                    times = ds["carbon_time"].var[:]
                    @test times == [0.0, 86400.0, 90000.0]
                    @test ds["carbon_iteration"][:] == [0, 96, 100]
                    uptake = ds["cumulative_air_sea_carbon_uptake_pgc"][:]
                    mean_rate = ds["mean_air_sea_carbon_uptake_rate_pgc_s"][:]
                    @test mean_rate[1] == 0
                    @test mean_rate[2:end] .* diff(times) ≈ diff(uptake)
                    if mode == :prescribed
                        @test all(iszero, ds["cumulative_carbon_emissions_pgc"][:])
                    else
                        @test all(iszero, ds["cumulative_implied_carbon_emissions_pgc"][:])
                    end
                end
            end
        end
    end
    @test length(variable_sets) == 2
    @test variable_sets[1] == variable_sets[2]
    @test length(first(variable_sets)) == 13
end
println("CARBON_DIAGNOSTICS_TEST_PASS")
end
