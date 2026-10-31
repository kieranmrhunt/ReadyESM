using ReadyESM
using Test
using SHA
using Dates

const R = ReadyESM
const O = R.Oceananigans

length(ARGS) in (3, 4) || error(
    "usage: seasonal_albedo_wiring_smoke.jl CONFIG ALBEDO_PATH ALBEDO_SHA256 [cpu|gpu|baseline_control]",
)
config_path, albedo_path, expected_sha = ARGS
mode = length(ARGS) == 4 ? Symbol(ARGS[4]) : :cpu
baseline_control = mode == :baseline_control
device = baseline_control ? :cpu : mode
device in (:cpu, :gpu) || baseline_control || error("device must be cpu, gpu or baseline_control")
isfile(albedo_path) || error("missing seasonal albedo file $albedo_path")
bytes2hex(sha256(read(albedo_path))) == expected_sha || error("seasonal albedo hash mismatch")

base = R.load_config(config_path)
settings = (; (name => getproperty(base, name) for name in fieldnames(typeof(base)))...)
config = R.ExperimentConfig(; merge(settings, (
    name = "seasonal_albedo_wiring_$(device)_smoke",
    device,
    start_date = DateTime(1993, 7, 1),
    duration_days = 0.25,
    output_dir = joinpath(pwd(), "seasonal_albedo_cpu_smoke_output"),
    terrarium_surface_albedo = :era5_spectral,
    terrarium_albedo_cycle_path = abspath(albedo_path),
    terrarium_albedo_cycle_sha256 = expected_sha,
    terrarium_albedo_crossing_band = :equal_halves,
))...)

simulation = R.build_dynamic_esm(config; steps = 1)
O.run!(simulation)

atmosphere = simulation.model.atmosphere
model = atmosphere.model
vars = atmosphere.variables
radiation = model.longwave_radiation
controller = model.callbacks[:seasonal_albedo]
solver_direct = R.RRTMGP.direct_sw_surface_albedo(radiation.solver)
solver_diffuse = R.RRTMGP.diffuse_sw_surface_albedo(radiation.solver)

direct_bands_match = Array(solver_direct) == controller.direct
diffuse_bands_match = Array(solver_diffuse) == controller.diffuse
if baseline_control
    @test !direct_bands_match
    @test !diffuse_bands_match
    println(
        "SEASONAL_ALBEDO_UNPATCHED_CONTROL_EXPECTED_FAIL " *
        "solver_direct_bands_match=false solver_diffuse_bands_match=false",
    )
    exit()
end
@test direct_bands_match
@test diffuse_bands_match

fraction = vec(Array(model.land_sea_mask.mask.data))
land_columns = controller.indices
@test !isempty(land_columns)
@test maximum(abs.(controller.direct[:, land_columns] .- 0.4)) > 0.05

# At the same atmospheric and surface state, compare the solved TOA shortwave
# with the old scalar 0.4 land albedo, then solve again through the fixed path.
held = vars.parameterizations.land.albedo.data
fill!(held, eltype(held)(0.4))
R._update_rrtmgp_state!(radiation, vars, model)
R._update_fluxes!(radiation)
scalar_toa_sw = vec(Array(R.RRTMGP.sw_flux_up(radiation.solver))[end, :])
scalar_surface_up = vec(Array(R.RRTMGP.sw_flux_up(radiation.solver))[1, :])

R._solve_rrtmgp!(radiation, vars, model)
seasonal_toa_sw = vec(Array(R.RRTMGP.sw_flux_up(radiation.solver))[end, :])
seasonal_surface_up = vec(Array(R.RRTMGP.sw_flux_up(radiation.solver))[1, :])
seasonal_surface_down = vec(Array(R.RRTMGP.sw_flux_dn(radiation.solver))[1, :])
@test Array(solver_direct) == controller.direct
@test Array(solver_diffuse) == controller.diffuse
@test maximum(abs.(Array(held) .- 0.4)) > 0.05

# On a sunlit, nearly all-land column, the change in reflected surface SW
# should equal land_fraction × Δalbedo × incident surface SW, within the
# adapter's Float32 flux-closure tolerance.
held_land = vec(Array(held))
ocean_albedo = vec(Array(vars.parameterizations.ocean.albedo.data))
eligible = [i for i in land_columns if fraction[i] > 0.9 && seasonal_surface_down[i] > 10]
@test !isempty(eligible)
column = eligible[argmin([
    sum(controller.direct[:, i]) / size(controller.direct, 1) for i in eligible
])]
land_fraction = fraction[column]
seasonal_land_albedo = held_land[column]
scalar_mixed_albedo = (1 - land_fraction) * ocean_albedo[column] + land_fraction * 0.4
observed_surface_reflection_change = seasonal_surface_up[column] -
    scalar_mixed_albedo * seasonal_surface_down[column]
expected_surface_reflection_change = land_fraction * (seasonal_land_albedo - 0.4) *
    seasonal_surface_down[column]
@test seasonal_land_albedo < 0.35
@test expected_surface_reflection_change < 0
@test isapprox(
    observed_surface_reflection_change,
    expected_surface_reflection_change;
    rtol = 1.0e-4,
    atol = 0.02,
)

weights = Float64.(Array(R._global_point_weights(model.spectral_grid)))
land_weights = weights .* fraction
land_mean_toa_change = sum(land_weights .* (seasonal_toa_sw .- scalar_toa_sw)) /
                       sum(land_weights)
land_mean_surface_reflection_change =
    sum(land_weights .* (seasonal_surface_up .- scalar_surface_up)) /
    sum(land_weights)
toa_transmission_ratio = land_mean_toa_change / land_mean_surface_reflection_change
@test abs(land_mean_toa_change) > 1.0
@test land_mean_surface_reflection_change < 0
@test isapprox(
    land_mean_toa_change,
    0.85 * land_mean_surface_reflection_change;
    rtol = 0.30,
)

result_tag = device == :cpu ?
    "SEASONAL_ALBEDO_CPU_WIRING_PASS " : "SEASONAL_ALBEDO_GPU_WIRING_PASS "
println(
    result_tag *
    "solver_bands_match=true " *
    "held_land_albedo_changed=true " *
    "surface_reflection_closes=true " *
    "land_mean_surface_reflection_change_W_m2=$(land_mean_surface_reflection_change) " *
    "land_mean_toa_sw_change_W_m2=$(land_mean_toa_change) " *
    "implied_upward_transmission=$(toa_transmission_ratio)",
)
