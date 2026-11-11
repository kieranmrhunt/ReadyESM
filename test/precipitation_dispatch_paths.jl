module PrecipitationDispatchPathTests
using ReadyESM, Test
const SW = ReadyESM.SpeedyWeather

struct ControlledConvection <: SW.AbstractConvection end
struct ControlledCondensation <: SW.AbstractCondensation
    target::Float32
end

function SW.parameterization!(ij, vars, ::ControlledConvection, model)
    vars.tendencies.grid.humidity[ij, 1] -= 2f-8
    vars.tendencies.grid.humidity[ij, 2] += 1f-8
    vars.parameterizations.rain_rate[ij] += 1f-7
    vars.parameterizations.rain_rate_convection[ij] = 1f-7
    vars.parameterizations.rain_convection[ij] += model.time_stepping.Δt_sec * 1f-7
    return nothing
end

function SW.parameterization!(ij, vars, process::ControlledCondensation, model)
    # Equal layer tendencies give a column-average loss of target/10.
    vars.tendencies.grid.humidity[ij, 1] -= process.target / 10f0
    vars.tendencies.grid.humidity[ij, 2] -= process.target / 10f0
    vars.parameterizations.rain_rate_large_scale[ij] = 2f-7
    vars.parameterizations.snow_rate_large_scale[ij] = 1f-7
    vars.parameterizations.rain_rate[ij] += 2f-7
    vars.parameterizations.snow_rate[ij] += 1f-7
    vars.parameterizations.rain_large_scale[ij] += model.time_stepping.Δt_sec * 2f-7
    vars.parameterizations.snow_large_scale[ij] += model.time_stepping.Δt_sec * 1f-7
    return nothing
end

model = (
    geometry = (σ_levels_thick = Float32[0.5, 0.5],),
    planet = (gravity = 10f0,),
    atmosphere = (water_density = 1000f0,),
    time_stepping = (Δt_sec = 3600f0,),
)
function makevars(; rain=3.9f-6, convection=0f0, large_rain=0f0,
                  large_snow=0f0, snow=large_snow)
    zeros2 = () -> zeros(Float32, 1, 2)
    zeros1 = () -> zeros(Float32, 1)
    p = (
        rain_rate = Float32[rain], rain_rate_convection = Float32[convection],
        rain_rate_large_scale = Float32[large_rain],
        snow_rate_large_scale = Float32[large_snow], snow_rate = Float32[snow],
        rain_convection = zeros1(), rain_large_scale = zeros1(),
        snow_large_scale = zeros1(),
        convective_humidity_tendency = zeros2(),
        convective_temperature_tendency = zeros2(),
        convective_precipitation_mass_correction = zeros1(),
        convective_precipitation_surface_pressure = zeros1(),
        large_scale_condensation_humidity_tendency = zeros2(),
        large_scale_condensation_temperature_tendency = zeros2(),
        large_scale_precipitation_mass_correction = zeros1(),
        large_scale_precipitation_surface_pressure = zeros1(),
    )
    return (
        tendencies = (grid = (humidity = zeros2(), temperature = zeros2()),),
        parameterizations = p,
        grid = (pressure_prev = Float32[100000],),
    )
end

@testset "controlled precipitation source dispatch" begin
    conv = makevars()
    SW.parameterization!(1, conv, ReadyESM.NetColumnBettsMillerConvection{ControlledConvection}(ControlledConvection()), model)
    cp = conv.parameterizations
    net = cp.rain_rate_convection[1]
    gross = (3.9f-6 + 1f-7) - 3.9f-6 # stored Float32 increment
    @test net > 0f0
    @test net != gross
    @test cp.rain_rate[1] == Float32(3.9f-6 + net)
    @test cp.convective_precipitation_mass_correction[1] == net - gross
    @test cp.rain_convection[1] == 3600f0 * 1f-7 + 3600f0 * (net - gross)
    @test cp.convective_humidity_tendency[1, :] == Float32[-2f-8, 1f-8]
    println("CONTROLLED_CONVECTION_PASS gross=", gross, " net=", net)

    # Both rain and snow corrections are nonzero when target < gross snow.
    imp = makevars(rain=1f-7, convection=1f-7, large_rain=0f0,
                   large_snow=0f0, snow=0f0)
    target = 5f-8
    SW.parameterization!(1, imp, ReadyESM.NetColumnImplicitCondensation{ControlledCondensation}(ControlledCondensation(target)), model)
    ip = imp.parameterizations
    @test ip.rain_rate_large_scale[1] == 0f0
    @test isapprox(ip.snow_rate_large_scale[1], target; rtol=1e-6)
    @test ip.rain_rate[1] == ip.rain_rate_convection[1] + ip.rain_rate_large_scale[1]
    @test ip.snow_rate[1] == ip.snow_rate_large_scale[1]
    @test ip.large_scale_precipitation_mass_correction[1] == (0f0 - 2f-7) + (ip.snow_rate_large_scale[1] - 1f-7)
    @test ip.rain_large_scale[1] == 3600f0 * 2f-7 + 3600f0 * (0f0 - 2f-7)
    @test ip.snow_large_scale[1] == 3600f0 * 1f-7 + 3600f0 * (ip.snow_rate_large_scale[1] - 1f-7)
    println("CONTROLLED_IMPLICIT_BOTH_CORRECTIONS_PASS rain=", ip.rain_rate_large_scale[1], " snow=", ip.snow_rate_large_scale[1])

    # A second case proves the stored sum when corrected rain remains positive.
    imp2 = makevars(rain=1f-7, convection=1f-7, large_rain=0f0,
                    large_snow=0f0, snow=0f0)
    target2 = 4f-7
    SW.parameterization!(1, imp2, ReadyESM.NetColumnImplicitCondensation{ControlledCondensation}(ControlledCondensation(target2)), model)
    ip2 = imp2.parameterizations
    @test ip2.rain_rate_large_scale[1] > 0f0
    @test isapprox(ip2.rain_rate_large_scale[1] + ip2.snow_rate_large_scale[1], target2; rtol=1e-6)
    @test ip2.rain_rate[1] == ip2.rain_rate_convection[1] + ip2.rain_rate_large_scale[1]
    @test ip2.snow_rate[1] == ip2.snow_rate_large_scale[1]
    println("CONTROLLED_IMPLICIT_POSITIVE_RAIN_PASS rain=", ip2.rain_rate_large_scale[1])
end
@testset "large-scale phase clears stale convection before the next step" begin
    # SpeedyWeather resets only the combined rain rate at the start of a step.
    # The large-scale process runs before convection, so this component still
    # contains the preceding step's value when large-scale condensation begins.
    ordered = makevars(rain=0f0, convection=9f-7, large_rain=0f0,
                       large_snow=0f0, snow=0f0)
    p = ordered.parameterizations
    large = ReadyESM.NetColumnImplicitCondensation{ControlledCondensation}(
        ControlledCondensation(4f-7))
    conv = ReadyESM.NetColumnBettsMillerConvection{ControlledConvection}(
        ControlledConvection())
    for step in 1:2
        p.rain_rate[1] = 0f0
        ordered.tendencies.grid.humidity .= 0f0
        ordered.tendencies.grid.temperature .= 0f0
        SW.parameterization!(1, ordered, large, model)
        @test p.rain_rate_convection[1] == 0f0
        @test p.rain_rate[1] == p.rain_rate_large_scale[1]
        SW.parameterization!(1, ordered, conv, model)
        @test p.rain_rate_convection[1] > 0f0
        @test p.rain_rate[1] == p.rain_rate_convection[1] +
                                 p.rain_rate_large_scale[1]
    end
    println("ORDERED_PRECIPITATION_TWO_STEP_PASS")
end

println("CONTROLLED_PRECIPITATION_DISPATCH_PASS")

end # module PrecipitationDispatchPathTests
