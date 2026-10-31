using ReadyESM
using Test
using SHA
using Dates

const R = ReadyESM
const O = R.Oceananigans

# CPU smoke of the end-of-run diagnostics with the ERA5 spectral land albedo on and off. The V3c month failed in
# `collect_dynamic_diagnostics` because `_atmosphere_radiative_surface_albedo` required band-invariant albedos. This
# runs 18 coupled steps per arm, then collects, saves and validates the diagnostics, and checks the albedo field and its
# source attribute in both arms.
length(ARGS) == 4 || error("usage: albedo_diagnostic_smoke.jl CONFIG ALBEDO_PATH ALBEDO_SHA256 OUTPUT_ROOT")
config_path, albedo_path, expected_sha, output_root = ARGS
isfile(albedo_path) || error("missing seasonal albedo file $albedo_path")
bytes2hex(sha256(read(albedo_path))) == expected_sha || error("seasonal albedo hash mismatch")
ispath(output_root) && error("refusing to reuse $output_root")

base = R.load_config(config_path)
settings = (; (name => getproperty(base, name) for name in fieldnames(typeof(base)))...)

function arm_config(arm)
    common = (name = "albedo_diagnostic_$(arm)_smoke", device = :cpu, start_date = DateTime(1993, 7, 1),
        duration_days = 0.25, output_dir = joinpath(abspath(output_root), String(arm)))
    arm == :spectral || return R.ExperimentConfig(; merge(settings, common)...)
    return R.ExperimentConfig(; merge(settings, common, (
        terrarium_surface_albedo = :era5_spectral,
        terrarium_albedo_cycle_path = abspath(albedo_path),
        terrarium_albedo_cycle_sha256 = expected_sha,
        terrarium_albedo_crossing_band = :equal_halves,
    ))...)
end

for arm in (:constant, :spectral)
    config = arm_config(arm)
    simulation = R.build_dynamic_esm(config; steps = 18)
    O.run!(simulation)
    model = simulation.model.atmosphere.model
    radiation = model.longwave_radiation
    direct = Float64.(Array(R.RRTMGP.direct_sw_surface_albedo(radiation.solver)))
    diffuse = Float64.(Array(R.RRTMGP.diffuse_sw_surface_albedo(radiation.solver)))
    diagnostics = R.collect_dynamic_diagnostics(simulation, config)
    paths = R.save_dynamic_diagnostics(diagnostics, config; render_figure = false)
    R.validate_dynamic_diagnostics(diagnostics)
    albedo = diagnostics.atmosphere_radiative_surface_albedo
    source = diagnostics.metadata.atmosphere_radiative_surface_albedo_source
    @test isfile(paths.netcdf_path)
    @test all((0 .<= albedo) .& (albedo .<= 1))
    if arm == :constant
        # Default path: unchanged definition and attribute, band 1 exactly.
        @test source == "retained_rrtmgp_direct_diffuse_band_invariant_solver_state"
        @test maximum(abs, direct .- diffuse) <= 8eps(Float32)
        @test albedo == vec(direct[1, :])
    else
        land = model.callbacks[:seasonal_albedo].indices
        @test source == "retained_rrtmgp_unweighted_band_and_direct_diffuse_mean_solver_state"
        @test maximum(abs, direct .- diffuse) > 8eps(Float32)
        @test albedo == vec(sum(direct .+ diffuse; dims = 1)) ./ (2 * size(direct, 1))
        @test !isempty(land) && maximum(abs.(albedo[land] .- 0.4)) > 0.05
    end
    println("ALBEDO_DIAGNOSTIC_SMOKE_ARM_PASS arm=$arm source=$source columns=$(length(albedo)) " *
        "mean=$(sum(albedo) / length(albedo)) netcdf=$(paths.netcdf_path)")
end
println("ALBEDO_DIAGNOSTIC_SMOKE_PASS")
