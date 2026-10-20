# BS08 runoff enthalpy: focused package regressions. These are reduced copies of the frozen external qualification
# gates (bs08_stage_b_v1/tests); those gates remain the qualification record.
using Test, Random
const RB = ReadyESM
const OB = RB.Oceananigans

hull_ok(T, sources; tol = 1e-9) = minimum(sources) - tol <= T <= maximum(sources) + tol
bs08_zstar_grid() = OB.TripolarGrid(OB.CPU(), Float64; size = (16, 12, 2), halo = (3, 3, 3),
    z = OB.MutableVerticalDiscretization([-100.0, -50.0, 0.0]))
bs08_plain_grid() = OB.TripolarGrid(OB.CPU(), Float64; size = (16, 12, 2), halo = (3, 3, 3), z = (-100.0, 0.0))

"Land sources routed to a tripolar grid with the native router; the exchange built by the BS08 builder."
function bs08_routed_land(grid, nsrc; owns = true)
    wet = RB._wet_ocean_cells(grid)
    picks = round.(Int, range(1, length(wet.wet_i); length = nsrc))
    longitude = [wet.λc[wet.wet_i[n], wet.wet_j[n]] for n in picks]
    latitude = [wet.φc[wet.wet_i[n], wet.wet_j[n]] for n in picks]
    area = [1e10 * (1 + 0.37 * k) for k in 1:nsrc]
    routing = RB._build_runoff_routing(longitude, latitude, area, grid; receivers_per_source = 4)
    xy = OB.Field{OB.Center, OB.Center, Nothing}(OB.RectilinearGrid(OB.CPU(); size = (nsrc, 1), x = (0, 1), y = (0, 1),
                                                                     topology = (OB.Bounded, OB.Bounded, OB.Flat)))
    state = owns ? (runoff_pending_depth = deepcopy(xy), runoff_pending_enthalpy = deepcopy(xy),
                    surface_excess_enthalpy = deepcopy(xy)) : (runoff_pending_depth = deepcopy(xy),)
    pending_depth = OB.interior(state.runoff_pending_depth)
    surface_runoff = similar(pending_depth); fill!(surface_runoff, 0)
    exchange = RB._runoff_exchange_storage(state, pending_depth, surface_runoff, (model = (grid = grid,),))
    land = RB.TerrariumRunoffLand(surface_runoff, routing.contribution_column, routing.contribution_weight,
        routing.target_i, routing.target_j, routing.offsets, nsrc, length(routing.target_i),
        routing.receivers_per_source, routing.routing_weighting, sum(area), nothing, exchange)
    return land, state, area
end

function bs08_row(iteration; source_water = 1e3, receiver_water = source_water, source_h = 2e10, receiver_h = source_h,
                  released = 0.0, received = 0.0, inconsistent = 0)
    return (; iteration, time_seconds = 900.0 * iteration, source_water_kgs = source_water,
            receiver_water_kgs = receiver_water, source_enthalpy_w = source_h, receiver_enthalpy_w = receiver_h,
            native_sst_heat_w = 1e10, bs08_minus_native_heat_w = receiver_h - 1e10, zero_flow = iszero(source_water),
            cumulative_released_enthalpy_j = released, cumulative_received_enthalpy_j = received,
            surface_excess_enthalpy_j = 1e12, pending_runoff_enthalpy_j = 2e12, cumulative_clamp_enthalpy_j = 0.0,
            partition_inconsistent_columns = inconsistent)
end

@testset "BS08 runoff enthalpy" begin
    rng = MersenneTwister(808)
    @testset "scalar laws" begin
        cv = RB.BS08_SURFACE_CV
        taur = 3600.0
        for _ in 1:2_000
            S = rand(rng) < 0.4 ? 0.0 : 10.0^(-8 + 6rand(rng)); Ts = -5 + 30rand(rng); E = S * cv * Ts
            rainfall = rand(rng) < 0.2 ? 0.0 : 10.0^(-9 + 4rand(rng)); snowfall = rand(rng) < 0.5 ? 0.0 : rand(rng) * rainfall
            Tair = -20 + 45rand(rng)
            Hin, _ = RB.bs08_incoming_enthalpy_rate(rainfall, snowfall, Tair)
            k = 10.0^(-9 + 5rand(rng)); D = S > 0 ? S / taur : 0.0; I = S > 0 ? min(D, k) : min(rainfall, k)
            R = rainfall + D - I
            p = RB.bs08_partition_rates(S, E, rainfall, Hin, D, I, R)
            @test p.consistent
            @test abs((p.infiltration + p.runoff) - (Hin + p.drainage)) <= 1e-12 * (abs(Hin) + abs(p.drainage)) + 1e-300
            src = Float64[]; rainfall > snowfall && push!(src, max(Tair, 0.0)); snowfall > 0 && push!(src, 0.0)
            S > 0 && push!(src, Ts)
            I > 0 && @test hull_ok(p.infiltration / (cv * I), src)
            R > 0 && @test hull_ok(p.runoff / (cv * R), src)
        end
        @test !RB.bs08_partition_rates(0.01, 0.01 * cv * 5, 1e-6, 1e-6 * cv * 10, 0.01 / 3600, 1e-7, 1.0).consistent
        @test RB.bs08_incoming_enthalpy_rate(1e-6, 1e-6, 15.0)[1] == 0.0
        @test RB.bs08_incoming_enthalpy_rate(2e-6, 1e-6, 4.0)[1] ≈ 1e-6 * cv * 4.0
        cvl, cvi, rhoL = 4.2e6, 1.9e6, 1000.0 * 333550.0
        @test RB.bs08_phase_aware_enthalpy(0.01, 5.0, 1.0, cvl, cvi, rhoL) ≈ 0.01 * cvl * 5.0
        @test RB.bs08_phase_aware_enthalpy(0.01, -3.0, 0.0, cvl, cvi, rhoL) ≈ 0.01 * (cvi * -3.0 - rhoL)
        d, c = RB.bs08_exfiltration_enthalpy(0.004, 6.0, 0.25, cvl)
        @test d == c && d ≈ 0.001 * cvl * 6.0
        for _ in 1:2_000
            Pw = 10.0^(-9 + 6rand(rng)); Ph = Pw * cv * (-2 + 25rand(rng)); dw = rand(rng) * Pw
            a, b = RB.bs08_consume(Pw, Ph, dw)
            @test abs(a + b - Ph) <= 4eps(abs(Ph)) + 1e-300
            dw > 0 && @test isapprox(a / dw, Ph / Pw; rtol = 1e-12)
        end
        rho, cpo, J, Toc = 1026.0, 3991.0, 1e-4, 12.0
        @test isapprox(RB.bs08_land_heat_content_correction(J * cpo * Toc, J, Toc, rho, cpo), 0.0; atol = 1e-15)
        @test RB.bs08_land_heat_content_correction(J * cpo * (Toc + 5), J, Toc, rho, cpo) > 0
    end

    @testset "option-off dispatch and storage" begin
        off, _, _ = bs08_routed_land(bs08_zstar_grid(), 3; owns = false)
        @test keys(off.runoff_exchange) == (:pending_depth,)
        @test !(off isa RB._BS08RunoffLand)
        @test OB.prognostic_state(off).schema == "readiesm_native_runoff_exchange_v1"
        @test isnothing(RB._bs08_diagnostic_history(off))
        on, _, _ = bs08_routed_land(bs08_zstar_grid(), 3)
        @test on isa RB._BS08RunoffLand
        @test on.runoff_exchange.guard.state === :idle
        @test_throws ErrorException bs08_routed_land(bs08_plain_grid(), 3)
    end

    @testset "routing conserves water and enthalpy; guard applies once" begin
        grid = bs08_zstar_grid()
        land, _, area = bs08_routed_land(grid, 6)
        ex = land.runoff_exchange
        rates = [1e-7, 4e-7, 0.0, 2.5e-8, 3e-7, 1e-9]
        land.surface_runoff .= reshape(rates, :, 1, 1)
        ex.held_enthalpy_rate .= reshape(RB.BS08_SURFACE_CV .* 11.25 .* rates, :, 1, 1)
        J = OB.Field{OB.Center, OB.Center, Nothing}(grid)
        RB._scatter_terrarium_runoff_state!(J, land)
        RB._bs08_scatter_runoff_enthalpy!(land)
        source = sum(area .* RB.BS08_SURFACE_CV .* 11.25 .* rates)
        @test RB._budget_scalar(RB._budget_integral(ex.receiver_enthalpy_flux)) ≈ source rtol = 2e-6
        T = OB.CenterField(grid); OB.set!(T, (x, y, z) -> 2.0 + 20 * cosd(y))
        Tsurf = view(T.data, :, :, 2:2)
        Jh = OB.Field{OB.Center, OB.Center, Nothing}(grid)
        clock = (iteration = 3, time = 2700.0)
        @test_throws ErrorException RB._bs08_apply_receiver_enthalpy!(Jh.data, land, J.data, Tsurf, grid, 1026.0, 3991.0, clock)
        ex.guard.state = :assembled
        RB._bs08_apply_receiver_enthalpy!(Jh.data, land, J.data, Tsurf, grid, 1026.0, 3991.0, clock)
        kept = copy(Jh.data)
        @test_throws ErrorException RB._bs08_apply_receiver_enthalpy!(Jh.data, land, J.data, Tsurf, grid, 1026.0, 3991.0, clock)
        @test isequal(Jh.data, kept) && ex.guard.applications == 1
    end

    @testset "closure, zero flow and output validation" begin
        @test RB._bs08_relative_closure(0.0, 0.0) == 0.0
        @test isnan(RB._bs08_relative_closure(0.0, 1.0)) && isnan(RB._bs08_relative_closure(1.0, 0.0))
        @test RB._bs08_relative_closure(2.0, 2.0) == 0.0
        zero_row = bs08_row(0; source_water = 0.0, source_h = 0.0)
        @test isnan(RB._bs08_routed_temperature(zero_row))
        warm = bs08_row(4; source_water = 1e3, source_h = 1e3 * RB.BS08_SURFACE_CV / 1000 * 7.5)
        @test RB._bs08_routed_temperature(warm) ≈ 7.5
        maps = (source_area_m2 = [1.0, 2.0], surface_excess_enthalpy_jm2 = [0.0, 1.0],
                pending_runoff_enthalpy_jm2 = [0.0, 0.0], cumulative_clamp_enthalpy_jm2 = [0.0, 0.0],
                receiver_enthalpy_flux_wm2 = zeros(3, 2))
        good = (; bs08_diagnostics = (; schema = RB._BS08_HISTORY_SCHEMA, interval_seconds = 3600.0,
                rows = [zero_row, warm], in_transit_enthalpy_j = 0.0, maps...))
        @test isnothing(RB._validate_bs08_diagnostics(good))
        for bad_row in (bs08_row(4; receiver_h = 0.0), bs08_row(4; receiver_water = 0.0), bs08_row(4; inconsistent = 2),
                        merge(warm, (; zero_flow = true)), bs08_row(4; receiver_h = 4e10), bs08_row(4; receiver_water = 2e3),
                        merge(warm, (; pending_runoff_enthalpy_j = NaN)), merge(warm, (; time_seconds = 0.0)))
            bad = (; bs08_diagnostics = merge(good.bs08_diagnostics, (; rows = [zero_row, bad_row])))
            @test_throws ErrorException RB._validate_bs08_diagnostics(bad)
        end
        @test isnothing(RB._validate_bs08_diagnostics((; other = 1)))     # option off: nothing to validate
        mktempdir() do directory
            path = joinpath(directory, "bs08.nc")
            RB.NCDatasets.NCDataset(path, "c") do ds
                RB.NCDatasets.defDim(ds, "terrarium_land_column", 2)
                RB.NCDatasets.defDim(ds, "ocean_x", 3); RB.NCDatasets.defDim(ds, "ocean_y", 2)
                RB._write_bs08_diagnostics!(ds, good, ("ocean_x", "ocean_y"), "1993-01-01 00:00:00")
            end
            RB.NCDatasets.NCDataset(path) do ds
                @test ds["bs08_iteration"][:] == [0, 4]
                t = ds["bs08_routed_water_mean_temperature_c"][:]
                @test ismissing(t[1]) || isnan(t[1])
                @test t[2] ≈ 7.5
                @test ds["bs08_zero_flow"][:] == Int8[1, 0]
                @test ds.attrib["bs08_enthalpy_reference"] == "liquid water at 0 degC"
            end
        end
    end

    @testset "history checkpoint v3: validate before restore, continue past the endpoint" begin
        land, _, area = bs08_routed_land(bs08_zstar_grid(), 6)
        ex = land.runoff_exchange
        history = RB.BS08DiagnosticHistory(area, 4, -1, NaN, 0.0, 0.0, 0.0, NamedTuple[])
        RB._BS08_DIAGNOSTIC_HISTORIES[ex.guard] = history
        ex.pending_depth .= reshape([0.0, 1e-9, 2e-3, 0.0, 4e-4, 1e-6], :, 1, 1)
        ex.pending_enthalpy .= RB.BS08_SURFACE_CV .* 5.0 .* ex.pending_depth
        append!(history.rows, [bs08_row(0; source_water = 0.0, source_h = 0.0), bs08_row(4; released = 9.0, received = 8.0),
                               bs08_row(8; released = 19.0, received = 17.0)])
        history.last_iteration, history.last_time_seconds = 8, 7200.0
        history.in_effect_receiver_enthalpy_w = 2e10
        history.cumulative_released_enthalpy_j, history.cumulative_received_enthalpy_j = 19.0, 17.0
        saved = OB.prognostic_state(land)
        @test saved.schema == "readiesm_native_runoff_exchange_v3"
        @test saved.diagnostic_history.rows == history.rows
        # an invalid history is refused before any exchange or history change
        snapshot = deepcopy(OB.prognostic_state(land))
        broken = merge(saved, (; diagnostic_history = merge(saved.diagnostic_history,
            (; rows = [bs08_row(4), bs08_row(8)]))))                     # missing launch row
        @test_throws ErrorException OB.restore_prognostic_state!(land, broken)
        @test OB.prognostic_state(land) == snapshot
        wrong_cadence = merge(saved, (; diagnostic_history = merge(saved.diagnostic_history, (; sample_every_n_steps = 96))))
        @test_throws ErrorException OB.restore_prognostic_state!(land, wrong_cadence)
        @test OB.prognostic_state(land) == snapshot
        # a valid restore keeps the history object and its full ledger
        empty!(history.rows); history.last_iteration = -1; history.cumulative_released_enthalpy_j = NaN
        fill!(ex.pending_enthalpy, NaN)
        OB.restore_prognostic_state!(land, saved)
        @test RB._bs08_diagnostic_history(land) === history
        @test OB.prognostic_state(land) == saved
        @test history.last_iteration == 8 && history.cumulative_released_enthalpy_j == 19.0
        push!(history.rows, bs08_row(12; released = 29.0, received = 27.0))   # continues past the saved endpoint
        @test length(history.rows) == 4 && history.rows[end].iteration == 12
        # invalid saved histories are refused before any exchange or history change (Codex 15:35 probe cases)
        populated = deepcopy(OB.prognostic_state(land))
        endpoint = saved.diagnostic_history.rows[end]
        for broken_history in (
            merge(saved.diagnostic_history, (; rows = [saved.diagnostic_history.rows[1:2]...,
                                                        merge(endpoint, (; pending_runoff_enthalpy_j = NaN))])),
            merge(saved.diagnostic_history, (; rows = [saved.diagnostic_history.rows[1:2]...,
                                                        merge(endpoint, (; time_seconds = -1.0))])),
            merge(saved.diagnostic_history, (; last_time_seconds = -1.0)),
            merge(saved.diagnostic_history, (; rows = [saved.diagnostic_history.rows[1:2]...,
                                                        merge(endpoint, (; cumulative_released_enthalpy_j = 99.0))])),
            merge(saved.diagnostic_history, (; in_effect_receiver_enthalpy_w = 1.0)),
            merge(saved.diagnostic_history, (; rows = [saved.diagnostic_history.rows[1:2]...,
                                                        merge(endpoint, (; receiver_enthalpy_w = 4e10))])),
        )
            @test_throws ErrorException OB.restore_prognostic_state!(land, merge(saved, (; diagnostic_history = broken_history)))
            @test OB.prognostic_state(land) == populated
        end
        # an invalid legacy restore leaves a populated destination unchanged
        v2 = (; schema = "readiesm_native_runoff_exchange_v2", runoff_exchange = saved.runoff_exchange,
              runoff_enthalpy = saved.runoff_enthalpy)
        bad_v2 = merge(v2, (; runoff_enthalpy = (pending_enthalpy_jm2 = fill(NaN, size(ex.pending_enthalpy)),
                                                 held_enthalpy_rate_wm2 = saved.runoff_enthalpy.held_enthalpy_rate_wm2)))
        @test_throws ErrorException OB.restore_prognostic_state!(land, bad_v2)
        @test OB.prognostic_state(land) == populated
        # a valid legacy v2 restore returns the history to launch state, so the launch-only rule applies
        OB.restore_prognostic_state!(land, v2)
        @test Array(ex.pending_enthalpy) == saved.runoff_enthalpy.pending_enthalpy_jm2
        @test isempty(history.rows) && history.last_iteration == -1 && history.cumulative_released_enthalpy_j == 0.0
        # the same holds for v1 and for a v3 checkpoint saved without history
        OB.restore_prognostic_state!(land, saved)
        @test history.last_iteration == 8
        OB.restore_prognostic_state!(land, (; schema = "readiesm_native_runoff_exchange_v1",
                                           runoff_exchange = saved.runoff_exchange))
        @test isempty(history.rows) && history.last_iteration == -1
        OB.restore_prognostic_state!(land, saved)
        OB.restore_prognostic_state!(land, merge(saved, (; diagnostic_history = nothing)))
        @test isempty(history.rows) && history.last_iteration == -1
    end
end
println("BS08_RUNOFF_ENTHALPY_PACKAGE_TESTS_PASS")
