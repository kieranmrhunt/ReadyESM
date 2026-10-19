@testset "Finite multi-plume configuration" begin
    R = ReadyESM
    S = R.SpeedyWeather
    grid = S.SpectralGrid(; NF=Float32, trunc=7, nlayers=8)
    base = R.ExperimentConfig()
    fields = (;(n=>getfield(base,n) for n in fieldnames(typeof(base)))...)
    options = (;atmosphere_convection=:finite_multiplume_c08_v1,
        atmosphere_saturation_law=:mixed_phase,
        atmosphere_convection_timescale_seconds=7200)
    config = R.ExperimentConfig(;merge(fields,options)...)
    @test R.validate(config;check_input_files=false) === config
    scheme = R._atmosphere_convection(config,grid)
    @test scheme isa R.FiniteMultiPlumeConvection
    @test scheme.spectrum.entrainment_m == Float32.(0.75 .* [0.3e-3,1e-3,3e-3])
    @test scheme.spectrum.shared.time_scale_s === 7200f0
    @test scheme.spectrum.shared.precip_conversion_m === Float32(0.0015)
    @test scheme.spectrum.shared.condensate_threshold === Float32(0.0005)
    @test scheme.saturation_law === :mixed_phase
    @test all(v->v isa S.ParameterizationVariable,S.variables(scheme))
    names = [v.name for v in S.variables(scheme)]
    @test length(names) == length(unique(names))
    @test :snow_rate_convection in names
    @test all(n->n in names,R.FINITE_MULTIPLUME_FIELDS)
    @test R._atmosphere_convection(base,grid) isa R.NetColumnBettsMillerConvection
    @test isempty(R._finite_multiplume_metadata(base))
    @test R._finite_multiplume_metadata(config).atmosphere_convection_solver_precision == "Float64"
    for field in (:atmosphere_c08_entrainment_scale,
                  :atmosphere_c08_precipitation_conversion_m_inv,
                  :atmosphere_c08_condensate_threshold_kg_kg)
        for value in (-1.0,Inf,NaN)
            bad = R.ExperimentConfig(;merge(fields,options,NamedTuple{(field,)}((value,)))...)
            @test_throws ArgumentError R.validate(bad;check_input_files=false)
        end
    end
    gpu = R.ExperimentConfig(;merge(fields,options,(;device=:gpu))...)
    @test_throws ArgumentError R.validate(gpu;check_input_files=false)
    prognostic = R.ForcingConfig(radiation=:rrtmgp_all_sky,
        cloud_scheme=:prognostic_condensate,cloud_condensate_retention_fraction=0.5)
    incompatible = R.ExperimentConfig(;merge(fields,options,(;forcing=prognostic))...)
    @test_throws ArgumentError R.validate(incompatible;check_input_files=false)
    grid64 = S.SpectralGrid(;NF=Float64,trunc=7,nlayers=8)
    @test_throws ArgumentError R.FiniteMultiPlumeConvection(grid64)
    for kw in ((;time_scale_seconds=Inf),(;time_scale_seconds=0),
               (;entrainment_scale=-1),(;precipitation_conversion_m_inv=NaN),
               (;condensate_threshold_kg_kg=Inf),(;saturation_law=:unknown))
        @test_throws ArgumentError R.FiniteMultiPlumeConvection(grid;kw...)
    end
    mktempdir() do path
        yaml = joinpath(path,"c08.yml")
        write(yaml,"""
        device: cpu
        atmosphere_convection: finite_multiplume_c08_v1
        atmosphere_convection_timescale_seconds: 3600
        atmosphere_c08_entrainment_scale: 0.9
        atmosphere_c08_precipitation_conversion_m_inv: 0.00075
        atmosphere_c08_condensate_threshold_kg_kg: 0.001
        """)
        loaded = R.load_config(yaml;check_input_files=false)
        configured = R._atmosphere_convection(loaded,grid)
        @test configured.spectrum.entrainment_m == Float32.(0.9 .* [0.3e-3,1e-3,3e-3])
        @test configured.spectrum.shared.time_scale_s === 3600f0
        @test configured.spectrum.shared.precip_conversion_m === Float32(0.00075)
        @test configured.spectrum.shared.condensate_threshold === Float32(0.001)
    end
end
println("FINITE_MULTIPLUME_CONFIGURATION_PASS")
