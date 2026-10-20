@testset "RRTMGP CH4 and N2O options" begin
    R = ReadyESM
    production_path = joinpath(R.PROJECT_ROOT, "config", "production.yml")
    production = load_config(production_path; check_input_files = false)
    @test isnothing(R.ForcingConfig().ch4_ppb) && isnothing(R.ForcingConfig().n2o_ppb)
    @test isnothing(production.forcing.ch4_ppb) && isnothing(production.forcing.n2o_ppb)

    raw = R.YAML.load_file(production_path)
    raw["forcing"]["ch4_ppb"] = 1746.2545166015625
    raw["forcing"]["n2o_ppb"] = 310.522216796875
    path = tempname() * ".yml"
    R.YAML.write_file(path, raw)
    pd = load_config(path; check_input_files = false)
    @test pd.forcing.ch4_ppb == 1746.2545166015625
    @test pd.forcing.n2o_ppb == 310.522216796875
    @test R.validate(pd; check_input_files = false) === pd

    fields(x) = (; (n => getfield(x, n) for n in fieldnames(typeof(x)))...)
    with_forcing(config; kwargs...) = R.ExperimentConfig(; merge(fields(config),
        (; forcing = R.ForcingConfig(; merge(fields(config.forcing), (; kwargs...))...)))...)
    for bad in ((; ch4_ppb = 0.0), (; n2o_ppb = -1.0), (; ch4_ppb = NaN), (; n2o_ppb = Inf),
                (; radiation = :speedy_simplified, ch4_ppb = 800.0),
                (; radiation = :speedy_simplified, n2o_ppb = 270.0))
        @test_throws ArgumentError R.validate(with_forcing(production; bad...); check_input_files = false)
    end

    grid = R.SpeedyWeather.SpectralGrid(; NF = Float32, trunc = 7, nlayers = 8)
    base = R.RRTMGPRadiation(grid)
    @test isnothing(base.ch4_vmr) && isnothing(base.n2o_vmr)
    @test R.RRTMGP.volume_mixing_ratio(base.solver, "ch4") == 1.9f-6
    @test R.RRTMGP.volume_mixing_ratio(base.solver, "n2o") == 0.34f-6
    over = R.RRTMGPRadiation(grid; ch4_ppb = pd.forcing.ch4_ppb, n2o_ppb = pd.forcing.n2o_ppb)
    @test over.ch4_vmr === Float32(pd.forcing.ch4_ppb * 1.0e-9)
    @test over.n2o_vmr === Float32(pd.forcing.n2o_ppb * 1.0e-9)
    @test R.RRTMGP.volume_mixing_ratio(over.solver, "ch4") == over.ch4_vmr
    @test R.RRTMGP.volume_mixing_ratio(over.solver, "n2o") == over.n2o_vmr
    @test R.RRTMGP.volume_mixing_ratio(over.solver, "co2") == R.RRTMGP.volume_mixing_ratio(base.solver, "co2")
    only_ch4 = R.RRTMGPRadiation(grid; ch4_ppb = 798.7992553710938)
    @test isnothing(only_ch4.n2o_vmr)
    @test R.RRTMGP.volume_mixing_ratio(only_ch4.solver, "n2o") == 0.34f-6
    @test_throws ArgumentError R.RRTMGPRadiation(grid; ch4_ppb = 0)
    # Same standard profile, lower CH4 and N2O: less longwave absorption, so more outgoing longwave at the top.
    toa_up(radiation) = sum(Array(R.RRTMGP.lw_flux_up(radiation.solver))[end, :])
    @test toa_up(over) > toa_up(base)
    @test toa_up(R.RRTMGPRadiation(grid; ch4_ppb = 1900.0, n2o_ppb = 340.0)) == toa_up(base)

    default_provenance = R._forcing_provenance(base, production)
    @test !any(k -> occursin("ch4", String(k)) || occursin("n2o", String(k)), keys(default_provenance))
    override_provenance = R._forcing_provenance(over, pd)
    @test collect(keys(override_provenance))[1:length(default_provenance)] == collect(keys(default_provenance))
    @test override_provenance.requested_ch4_ppb == pd.forcing.ch4_ppb
    @test override_provenance.applied_ch4_ppb == Float64(over.ch4_vmr) * 1.0e9
    @test override_provenance.applied_n2o_ppb == Float64(over.n2o_vmr) * 1.0e9
    @test override_provenance.ch4_forcing_mode == "prescribed_well_mixed_volume_mixing_ratio"
end
println("RRTMGP_TRACE_GASES_UNIT_PASS")
