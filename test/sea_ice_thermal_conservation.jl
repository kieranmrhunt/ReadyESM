module SeaIceThermalConservationTest

using ReadyESM, Test
const O = ReadyESM.Oceananigans
const C = ReadyESM.ClimaSeaIce
const IT = C.SeaIceThermodynamics

@testset "Zero-tendency ice mass flux" begin
    thermodynamics = (concentration_evolution=IT.ProportionalEvolution(),)
    for (FT, DT) in ((Float32, Float32), (Float32, Float64), (Float64, Float64)),
        h in FT.((0.12345, 1.2345, 4.321, 14.987)),
        a in FT.((0.0012345, 0.12345, 0.54321, 0.98765)),
        dt in DT.((300, 450, 900))
        hn, an = IT.ice_volume_update(thermodynamics, zero(FT), h, a, FT(0.05), dt)
        # Use the kernel's pre-storage arithmetic, not rounded h/A. Promoting
        # these helper returns can create a flux while the stored state agrees.
        flux = FT(FT(900) * (hn * an - h * a) / dt)
        @test iszero(flux)
        @test isequal(FT(hn), h)
        @test isequal(FT(an), a)
    end
end

@testset "Native layered thermodynamics without forcing" begin
    for (FT, DT) in ((Float32, Float32), (Float32, Float64), (Float64, Float64))
        grid = O.RectilinearGrid(O.CPU(), FT;
            size=(2, 2, 1), extent=(2, 2, 1),
            topology=(O.Periodic, O.Periodic, O.Bounded))
        top = IT.PrescribedTemperature(zero(FT))
        bottom = IT.IceWaterThermalEquilibrium(zero(FT))
        ice = IT.sea_ice_slab_thermodynamics(grid;
            top_heat_boundary_condition=top, bottom_heat_boundary_condition=bottom)
        snow = IT.snow_slab_thermodynamics(grid;
            top_heat_boundary_condition=top, bottom_heat_boundary_condition=bottom)
        model = C.SeaIceModel(grid; ice_thermodynamics=ice, snow_thermodynamics=snow,
            top_heat_flux=zero(FT), bottom_heat_flux=zero(FT), snowfall=zero(FT),
            dynamics=nothing, advection=nothing, timestepper=:ForwardEuler)
        # Binary-exact initial volumes; no flooding or external heat/water.
        O.set!(model; h=FT(1), ℵ=FT(0.5), hs=FT(0.125))
        IT.thermodynamic_time_step!(model, model.ice_thermodynamics,
            model.snow_thermodynamics, DT(900))
        @test all(==(FT(1)), O.interior(model.ice_thickness))
        @test all(==(FT(0.5)), O.interior(model.ice_concentration))
        @test all(==(FT(0.125)), O.interior(model.snow_thickness))
        @test all(iszero, O.interior(model.mass_fluxes.thermodynamics.ice))
        @test all(iszero, O.interior(model.mass_fluxes.thermodynamics.snow))
        @test all(iszero, O.interior(model.mass_fluxes.intercepted_snowfall))
        @test all(isfinite, O.interior(model.ice_thermodynamics.top_surface_temperature))
        @test all(isfinite(model.snow_thermodynamics.top_surface_temperature[i, j, 1])
                  for i in 1:2, j in 1:2)
    end
end

end
