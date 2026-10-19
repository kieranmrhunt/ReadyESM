@testset "Conservative convective momentum configuration" begin
    R = ReadyESM
    S = R.SpeedyWeather
    grid = S.SpectralGrid(; NF=Float32, trunc=7, nlayers=8)
    native = R.ExperimentConfig()
    fields = (;(n=>getfield(native,n) for n in fieldnames(typeof(native)))...)
    options = (;atmosphere_convection=:finite_multiplume_c07_v1,
        atmosphere_convection_timescale_seconds=7200)
    config = R.ExperimentConfig(;merge(fields,options)...)
    @test R.validate(config;check_input_files=false) === config
    scheme = R._atmosphere_convection(config,grid)
    @test scheme isa R.C07HostAdapter.C07ThermalReturn
    @test scheme.enabled
    @test scheme.parent.base isa R.FiniteMultiPlumeConvection
    @test scheme.parent.strength === 0.5f0
    @test scheme.parent.maximum_subcycles == 16
    names = [v.name for v in S.variables(scheme)]
    @test length(names) == length(unique(names))
    @test all(n->n in names,(:snow_rate_convection,:c07_temperature_tendency,
        :c07_u_tendency,:c07_v_tendency,:c07_total_energy_residual))
    @test R._atmosphere_convection(native,grid) isa R.NetColumnBettsMillerConvection
    @test isempty(R._convective_momentum_metadata(native))
    @test isempty(R._finite_multiplume_metadata(native))
    @test R._finite_multiplume_metadata(config).atmosphere_convection_solver_precision == "Float64"
    @test R._convective_momentum_metadata(config).atmosphere_c07_momentum_strength == 0.5
    for override in ((;device=:gpu), (;atmosphere_c07_momentum_strength=-1.),
                     (;atmosphere_c07_momentum_strength=Inf), (;atmosphere_c07_momentum_strength=NaN),
                     (;atmosphere_c07_momentum_strength=1e100),
                     (;atmosphere_c07_maximum_subcycles=0), (;atmosphere_c07_maximum_subcycles=-1))
        bad = R.ExperimentConfig(;merge(fields,options,override)...)
        @test_throws ArgumentError R.validate(bad;check_input_files=false)
    end
    prognostic = R.ForcingConfig(radiation=:rrtmgp_all_sky,
        cloud_scheme=:prognostic_condensate,cloud_condensate_retention_fraction=0.5)
    @test_throws ArgumentError R.validate(R.ExperimentConfig(;
        merge(fields,options,(;forcing=prognostic))...);check_input_files=false)
    for override in ((;strength=-1.),(;strength=Inf),(;strength=1e100),(;maximum_subcycles=0))
        @test_throws ArgumentError R.ConvectiveMomentumTransport(scheme.parent.base;override...)
    end
    zero = R.ConvectiveMomentumTransport(scheme.parent.base;strength=0)
    @test iszero(zero.parent.strength) && zero.enabled
    mktempdir() do path
        file = joinpath(path,"c07.yml")
        write(file,"""
        device: cpu
        atmosphere_convection: finite_multiplume_c07_v1
        atmosphere_convection_timescale_seconds: 3600
        atmosphere_c08_precipitation_conversion_m_inv: 0.00075
        atmosphere_c08_condensate_threshold_kg_kg: 0.001
        atmosphere_c07_momentum_strength: 0.25
        atmosphere_c07_maximum_subcycles: 8
        """)
        loaded = R.load_config(file;check_input_files=false)
        configured = R._atmosphere_convection(loaded,grid)
        @test configured.parent.strength === 0.25f0
        @test configured.parent.maximum_subcycles == 8
        @test configured.parent.base.spectrum.shared.time_scale_s === 3600f0
        @test configured.parent.base.spectrum.shared.precip_conversion_m === Float32(.00075)
        @test configured.parent.base.spectrum.shared.condensate_threshold === Float32(.001)
    end
end
println("C07_MOMENTUM_CONFIGURATION_PASS")
