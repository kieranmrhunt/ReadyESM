#=
Spectral truncation of the atmosphere's prognostic state.

SpeedyWeather keeps one degree more than the truncation in its spectral arrays (l = 0 ... lmax + 1) for vector
operations. Its own initial conditions and orography zero that extra row with spectral_truncation!
("set the lmax+1 harmonics to zero"), and the leapfrog truncates every TENDENCY, but never the prognostic variables
themselves. ReadyESM's ERA5 initial conditions transform ln ps, temperature and humidity from grid fields without
truncating, so the launch leaves non-zero l = lmax + 1 coefficients that are never touched again: a frozen grid-scale
pattern for the whole run. In D (P3/P3b restarts) ln ps carries 1258x more power at l = lmax + 1 than at l = lmax,
temperature 20x and humidity 43x, unchanged from day 120 to day 240. In a standalone T31/L27 model with D's humidity
path it produces a non-process atmospheric water source of +0.0147 kg m-2/day; truncating once after the launch
removes it (+0.0001). (blue_sky_20261008/water_source_v1.)

`truncate_prognostic_lmax_plus_1!` zeroes that row for every prognostic spectral field and leapfrog level. Coefficients
with l <= lmax are untouched, so states already truncated are unchanged bitwise.
=#

const _TRUNCATED_PROGNOSTICS = (:vorticity, :divergence, :temperature, :humidity, :pressure)

"""Zero the l = lmax + 1 row of every prognostic spectral field (all leapfrog levels, tracers included)."""
function truncate_prognostic_lmax_plus_1!(vars::SpeedyWeather.Variables)
    prognostic = vars.prognostic
    for name in _TRUNCATED_PROGNOSTICS
        hasproperty(prognostic, name) || continue
        variable = getproperty(prognostic, name)
        for lf in 1:2
            SpeedyWeather.SpeedyTransforms.spectral_truncation!(SpeedyWeather.get_step(variable, lf))
        end
    end
    if hasproperty(prognostic, :tracers)
        for variable in values(prognostic.tracers), lf in 1:2
            SpeedyWeather.SpeedyTransforms.spectral_truncation!(SpeedyWeather.get_step(variable, lf))
        end
    end
    return vars
end

"""Apply the launch truncation when `config.atmosphere_initial_spectral_truncation` is set, after the initial
conditions are in spectral space and before the first grid transform (SpeedyWeather simulation initialize)."""
function _initial_spectral_truncation!(simulation, config::ExperimentConfig)
    config.atmosphere_initial_spectral_truncation && truncate_prognostic_lmax_plus_1!(simulation.variables)
    return simulation
end

const _GRID_PREVIOUS = (:u_prev, :v_prev, :temperature_prev, :pressure_prev, :humidity_prev)

function _copy_grid_previous(vars::SpeedyWeather.Variables)
    grid = vars.grid
    saved = Dict{Symbol, Any}(name => copy(getproperty(grid, name)) for name in _GRID_PREVIOUS if hasproperty(grid, name))
    if hasproperty(grid, :tracers)
        for (name, field) in pairs(grid.tracers)
            endswith(String(name), "_prev") && (saved[Symbol(:tracer_, name)] = copy(field))
        end
    end
    return saved
end

function _restore_grid_previous!(vars::SpeedyWeather.Variables, saved)
    grid = vars.grid
    for name in _GRID_PREVIOUS
        haskey(saved, name) && (getproperty(grid, name) .= saved[name])
    end
    if hasproperty(grid, :tracers)
        for (name, field) in pairs(grid.tracers)
            key = Symbol(:tracer_, name)
            haskey(saved, key) && (field .= saved[key])
        end
    end
    return vars
end

"""
Remove the l = lmax + 1 row from a RESTORED coupled atmosphere (e.g. a checkpoint of a run launched without the launch
truncation) and recompute its current grid state from the truncated spectra. This changes the state: a continuation
after it is not an exact restart.

A restore copies the grid state, including the `*_prev` grids (the previous step's state, used by the
parameterizations), without any transform. The recompute here uses `transform!(...; initialize = true)`, which would also
overwrite `*_prev` with the current state. With `preserve_previous = true` (default) the restored `*_prev` grids are put
back bit for bit, so the row removal (and the current grid recomputed from it) is the only change. `preserve_previous =
false` keeps the initialisation semantics (`*_prev` = current).

SpeedyWeather's transform! derives u and v from radius-scaled vorticity and divergence, the state of a stepping
atmosphere (the coupled model calls `SpeedyWeather.timestep!` and never unscales between steps; `run!`'s finalize does).
An unscaled state would give u/radius, so it is refused.
"""
function truncate_restored_atmosphere_lmax_plus_1!(simulation; preserve_previous::Bool = true)
    atmosphere = simulation.model.atmosphere
    vars = atmosphere.variables
    radius = atmosphere.model.planet.radius
    vars.prognostic.scale[] == radius || throw(ArgumentError(
        "truncate_restored_atmosphere_lmax_plus_1! needs the radius-scaled state of a stepping atmosphere " *
        "(prognostic.scale[] == planet radius = $radius), got scale $(vars.prognostic.scale[])",
    ))
    saved = preserve_previous ? _copy_grid_previous(vars) : nothing
    truncate_prognostic_lmax_plus_1!(vars)
    SpeedyWeather.transform!(vars, 2, atmosphere.model; initialize = true)
    preserve_previous && _restore_grid_previous!(vars, saved)
    return simulation
end
