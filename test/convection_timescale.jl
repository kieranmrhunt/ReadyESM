using Dates

@testset "Convection timescale" begin
    R = ReadyESM
    default = R.ExperimentConfig()
    @test default.atmosphere_convection_timescale_seconds == 14400
    grid = R.SpeedyWeather.SpectralGrid(; NF=Float32, trunc=7, nlayers=8)
    fields = (; (n => getfield(default,n) for n in fieldnames(typeof(default)))...)
    for seconds in (7200,14400), scheme in (:betts_miller_constant_rh,:betts_miller_sigma_rh_v1)
        c = R.ExperimentConfig(; merge(fields,(;
            atmosphere_convection=scheme,atmosphere_convection_timescale_seconds=seconds))...)
        @test R.validate(c;check_input_files=false) === c
        wrapper = R._atmosphere_convection(c,grid)
        @test wrapper isa R.NetColumnBettsMillerConvection
        convection = scheme == :betts_miller_constant_rh ? wrapper.convection : wrapper.convection.convection
        expected = R.SpeedyWeather.BettsMillerConvection(grid;time_scale=Second(seconds))
        for n in fieldnames(typeof(expected))
            @test isequal(getfield(convection,n),getfield(expected,n))
        end
    end
    actual = R._atmosphere_convection(default,grid).convection
    expected = R.SpeedyWeather.BettsMillerConvection(grid)
    for n in fieldnames(typeof(expected))
        @test isequal(getfield(actual,n),getfield(expected,n))
    end
    for seconds in (0,-1)
        c = R.ExperimentConfig(;merge(fields,(;atmosphere_convection_timescale_seconds=seconds))...)
        @test_throws ArgumentError R.validate(c;check_input_files=false)
    end
end
println("CONVECTION_TIMESCALE_UNIT_PASS")
