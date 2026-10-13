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

"""
Remove the l = lmax + 1 row from a RESTORED coupled atmosphere (e.g. a checkpoint of a run launched without the launch
truncation) and recompute its grid state from the truncated spectra, as at initialisation (grid and its `*_prev`
copy both from the current level). This changes the state: a continuation after it is not an exact restart.
"""
function truncate_restored_atmosphere_lmax_plus_1!(simulation)
    atmosphere = simulation.model.atmosphere
    truncate_prognostic_lmax_plus_1!(atmosphere.variables)
    SpeedyWeather.transform!(atmosphere.variables, 2, atmosphere.model; initialize = true)
    return simulation
end
