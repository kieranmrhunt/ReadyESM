using Dates
using ReadyESM

const R = ReadyESM
const O = R.Oceananigans

# Coupled CPU smoke of the opt-in large-scale snowmelt cooling on the production profile: T7/L8, ERA5 1 January 1993
# start, small analytic ocean, 900 s step, six hours, with the option off and on. Each arm collects, saves and
# validates the end-of-run production diagnostics, the step that failed in the first July month (M1, 63795338_3)
# because the validator did not know the provenance text of the option. With the option on, the provenance must name
# the snowmelt cooling and the passive snowmelt output must exist and be non-negative. This is not a climate test.
length(ARGS) == 1 || error("usage: large_scale_melt_cooling_smoke.jl OUTPUT_DIR")
output = abspath(ARGS[1])
ispath(output) && error("refusing to reuse $output")
fields(x) = (; (name => getfield(x, name) for name in fieldnames(typeof(x)))...)
era5 = "/gws/ssde/j25a/ncas_climate/vol1/users/kieran/readiesm/data/era5_initial_state"
base = R.load_config(joinpath(R.PROJECT_ROOT, "config", "production.yml"); check_input_files = false)
base.atmosphere_large_scale_precipitation == :net_column_conservative || error("production precipitation changed")

for melt in (false, true)
    name = melt ? "on" : "off"
    config = R.ExperimentConfig(; merge(fields(base), (
        name = "large_scale_melt_cooling_smoke_$name", device = :cpu, truncation = 7, nlayers = 8,
        duration_days = 0.25, analysis_start_day = 0.0, start_date = DateTime(1993, 1, 1),
        atmosphere_timestep_at_t31_minutes = 3.75,
        atmosphere_initial_conditions = :era5_instantaneous, era5_month = 1,
        era5_pressure_levels_path = joinpath(era5, "era5_pressure_levels_stratosphere_19930101T0000Z.nc"),
        era5_single_levels_path = joinpath(era5, "era5_single_levels_19930101T0000Z.nc"),
        era5_land_state_path = joinpath(era5, "era5_land_state_19930101T0000Z.nc"),
        ocean_initial_conditions = :analytic, ecco_initial_conditions_directory = "", ocean_grid = :idealized,
        ocean_dynamics = false, sea_ice_dynamics = false, ocean_nlongitude = 8, ocean_nlatitude = 4, ocean_nlayers = 1,
        output_dir = joinpath(output, name), atmosphere_large_scale_melt_cooling = melt,
    ))...)
    R.validate(config)
    simulation = R.build_dynamic_esm(config)
    lsc = simulation.model.atmosphere.model.large_scale_condensation
    (lsc.condensation isa R.MeltCoolingImplicitCondensation) == melt || error("$name: wrong large-scale scheme")
    O.run!(simulation)
    simulation.model.clock.iteration == 24 || error("$name: expected 24 coupled steps")
    diagnostics = R.collect_dynamic_diagnostics(simulation, config)
    paths = R.save_dynamic_diagnostics(diagnostics, config; render_figure = false)
    R.validate_dynamic_diagnostics(diagnostics)
    scheme = diagnostics.metadata.atmosphere_large_scale_precipitation_scheme
    if melt
        scheme == "implicit condensation net column humidity drying with snowmelt cooling" ||
            error("on: provenance does not name the snowmelt cooling: $scheme")
        rate = Array(simulation.model.atmosphere.variables.parameterizations.large_scale_snowmelt_rate.data)
        all(x -> isfinite(x) && x >= 0, rate) || error("on: invalid snowmelt output")
    else
        scheme == "implicit condensation net column humidity drying" || error("off: provenance changed: $scheme")
    end
    println("MELT_SMOKE_ARM arm=$name steps=24 diagnostics_validated=true provenance=\"$scheme\" netcdf=$(paths.netcdf_path)")
    flush(stdout)
end
println("LARGE_SCALE_MELT_COOLING_SMOKE_PASS")
flush(stdout)
