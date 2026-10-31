using Test
using ReadyESM

# Oracle O-007 and E-004: the standard van Genuchten–Mualem conductivity against an independent BigFloat reference,
# the released S^(n/(n+1)) form kept unchanged, the published ratios between them, the full
# hydraulic_conductivity path through a SoilComposition, and the production default.
function mualem_conductivity_test()
    R = ReadyESM
    T = R.Terrarium
    reference(S, n) = (S = big(S); n = big(n); m = 1 - 1 / n; sqrt(S) * (1 - (1 - S^(1 / m))^m)^2)
    released(S, n) = (S = big(S); n = big(n); sqrt(S) * (1 - (1 - S^(n / (n + 1)))^((n - 1) / n))^2)
    for NF in (Float64, Float32), n in (2.0, 1.5, 3.0), S in (0.05, 0.1, 0.2, 0.5, 0.9, 0.99, 1.0)
        value = T.vangenuchten_mualem_relative_conductivity(NF(S), NF(n))
        @test value isa NF
        @test Float64(value) ≈ Float64(reference(NF(S), NF(n))) rtol = (NF == Float64 ? 1e-12 : 2e-5)
    end
    # The released law over the corrected one for n = 2, as Oracle reported from the frozen calls: 8.55, 87.35, 519.49.
    for (S, ratio) in ((0.5, 8.55), (0.2, 87.35), (0.1, 519.49))
        @test Float64(released(S, 2) / reference(S, 2)) ≈ ratio rtol = 2e-3
    end
    @test T.vangenuchten_mualem_relative_conductivity(1.0, 2.0) == 1.0           # saturated endpoint unchanged

    # Through the actual hydraulic_conductivity path: unfrozen soil (no ice impedance) and a frozen case.
    for (name, cond) in (("released", T.UnsatKVanGenuchten(Float64)), ("mualem", T.UnsatKVanGenuchtenMualem(Float64)))
        hyd = T.ConstantSoilHydraulics(Float64; swrc = T.VanGenuchten(α = 2.0, n = 2.0), unsat_hydraulic_cond = cond)
        for S in (0.1, 0.2, 0.5, 1.0)
            soil = T.SoilComposition(porosity = 0.4, saturation = S, liquid = 1.0)
            Ksat = T.saturated_hydraulic_conductivity(hyd, soil.solid.texture)
            K = T.hydraulic_conductivity(hyd, soil)
            expected = name == "mualem" ? reference(S, 2.0) : released(S, 2.0)
            @test K / Ksat ≈ Float64(expected) rtol = 1e-12
        end
        frozen = T.SoilComposition(porosity = 0.4, saturation = 0.6, liquid = 0.5)
        Ksat = T.saturated_hydraulic_conductivity(hyd, frozen.solid.texture)
        Sliq = T.water(frozen) / T.porosity(frozen)
        base = name == "mualem" ? reference(Sliq, 2.0) : released(Sliq, 2.0)
        @test T.hydraulic_conductivity(hyd, frozen) / Ksat ≈ Float64(base) * 10^(-7 * 0.5) rtol = 1e-10
    end

    # Production keeps the released law, and the option selects the corrected type.
    production = R.load_config(joinpath(R.PROJECT_ROOT, "config/production.yml"); check_input_files = false)
    @test production.terrarium_unsat_conductivity == :van_genuchten_released
    @test R._prognostic_soil_hydrology(Float64).hydraulic_properties.unsat_hydraulic_cond isa T.UnsatKVanGenuchten
    @test R._prognostic_soil_hydrology(Float64, :van_genuchten_mualem).hydraulic_properties.unsat_hydraulic_cond isa
        T.UnsatKVanGenuchtenMualem
    @test_throws ArgumentError R.validate(R.ExperimentConfig(; terrarium_unsat_conductivity = :wrong))
end

@testset "van Genuchten–Mualem unsaturated conductivity (O-007)" begin
    mualem_conductivity_test()
end
