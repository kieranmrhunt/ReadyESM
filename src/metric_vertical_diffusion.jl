"""
    MetricBulkRichardsonDiffusion(spectral_grid; kwargs...)

Dimensionally consistent, explicitly stable bulk-Richardson vertical mixing
for SpeedyWeather's sigma-coordinate atmosphere.

SpeedyWeather diagnoses a physical kinematic diffusivity `K` in m² s⁻¹, while
its vertical finite-volume operator differentiates in dimensionless sigma.  The
physical density-weighted flux-form operator

    (1 / ρ) ∂z (ρ K ∂z x)

becomes

    ∂σ (K (ρg / ps)² ∂σ x)

after hydrostatic conversion, with `ρg / ps = σg / (R_d Tᵥ)`.  This wrapper
retains the upstream Richardson-number and physical-`K` diagnosis, applies that
metric conversion, and uses adjacent sigma levels with no flux at the diagnosed
boundary-layer top and at the surface.  The latter makes every diffused column
conservative under the model's sigma-layer weights.  Diffusion is evaluated
from the previous leapfrog state so that the prognostic update is a one-level
forward step, and a single column-wide scale factor enforces the monotonic
explicit finite-volume stability bound without changing the relative physical
mixing profile or breaking face-flux cancellation.
"""
struct MetricBulkRichardsonDiffusion{D} <:
       SpeedyWeather.AbstractVerticalDiffusion
    diffusion::D
end

Adapt.@adapt_structure MetricBulkRichardsonDiffusion

const ATMOSPHERE_VERTICAL_DIFFUSION_PROVENANCE =
    "metric_bulk_richardson_physical_K_sigma_flux_previous_state_cfl_v2"
const ATMOSPHERE_VERTICAL_DIFFUSION_CONTROL_PROVENANCE =
    "speedyweather_0.21.1_bulk_richardson_zero_adjacent_transport_control"
const METRIC_VERTICAL_DIFFUSION_COURANT_LIMIT = 0.9

function MetricBulkRichardsonDiffusion(
    spectral_grid::SpeedyWeather.SpectralGrid;
    kwargs...,
)
    diffusion = SpeedyWeather.BulkRichardsonDiffusion(
        spectral_grid;
        kwargs...,
    )
    return MetricBulkRichardsonDiffusion{typeof(diffusion)}(diffusion)
end

_atmosphere_vertical_diffusion(spectral_grid) =
    MetricBulkRichardsonDiffusion(spectral_grid)

function _atmosphere_vertical_diffusion(config::ExperimentConfig, spectral_grid)
    config.atmosphere_vertical_diffusion == :metric_bulk_richardson &&
        return MetricBulkRichardsonDiffusion(spectral_grid)
    config.atmosphere_vertical_diffusion ==
        :speedyweather_0_21_1_zero_operator_control &&
        return SpeedyWeather.BulkRichardsonDiffusion(spectral_grid)
    error(
        "unsupported atmosphere vertical diffusion " *
        "$(config.atmosphere_vertical_diffusion)",
    )
end

_atmosphere_vertical_diffusion_provenance(
    ::MetricBulkRichardsonDiffusion,
) = ATMOSPHERE_VERTICAL_DIFFUSION_PROVENANCE

_atmosphere_vertical_diffusion_provenance(
    ::SpeedyWeather.BulkRichardsonDiffusion,
) = ATMOSPHERE_VERTICAL_DIFFUSION_CONTROL_PROVENANCE

SpeedyWeather.variables(diffusion::MetricBulkRichardsonDiffusion) = (
    SpeedyWeather.variables(diffusion.diffusion)...,
    SpeedyWeather.ParameterizationVariable(
        :vertical_diffusion_cfl_scale,
        SpeedyWeather.Grid2D();
        desc = "Column scale applied to preserve explicit vertical-diffusion stability",
        units = "1",
    ),
)

SpeedyWeather.initialize!(
    diffusion::MetricBulkRichardsonDiffusion,
    model::SpeedyWeather.AbstractModel,
) = SpeedyWeather.initialize!(diffusion.diffusion, model)

Base.@propagate_inbounds function _metric_vertical_diffusion!(
    ij,
    tendency,
    variable,
    sigma_diffusivity,
    boundary_layer_top,
    diffusion,
)
    (; ∇²_above, ∇²_below) = diffusion
    nlayers = size(tendency, 2)

    for k in boundary_layer_top:nlayers
        # The diagnosed K is identically zero above boundary_layer_top.  Making
        # its upper face impermeable avoids a one-sided flux into a layer whose
        # tendency is not updated.  The last layer similarly has no surface
        # flux.  Interior face contributions cancel exactly after multiplying
        # by sigma-layer thickness.
        k_above = k == boundary_layer_top ? k : k - 1
        k_below = k == nlayers ? k : k + 1
        flux_below =
            (variable[ij, k_below] - variable[ij, k]) *
            (sigma_diffusivity[ij, k_below] + sigma_diffusivity[ij, k])
        flux_above =
            (variable[ij, k] - variable[ij, k_above]) *
            (sigma_diffusivity[ij, k] + sigma_diffusivity[ij, k_above])
        tendency[ij, k] +=
            ∇²_below[k] * flux_below - ∇²_above[k] * flux_above
    end
    return nothing
end

"""
    _limit_metric_sigma_diffusivity!(ij, Kσ, k_top, diffusion, dt)

Scale all face diffusivities in one column by the same factor when necessary
to make the explicit finite-volume update monotonic over `dt` seconds.  A
uniform column scale preserves the symmetric interior face fluxes and hence
the sigma-mass-weighted conservation property of `_metric_vertical_diffusion!`.
"""
Base.@propagate_inbounds function _limit_metric_sigma_diffusivity!(
    ij,
    sigma_diffusivity,
    boundary_layer_top,
    diffusion,
    maximum_timestep_seconds,
)
    (; ∇²_above, ∇²_below) = diffusion
    nlayers = size(sigma_diffusivity, 2)
    boundary_layer_top > nlayers && return one(maximum_timestep_seconds)
    reference = sigma_diffusivity[ij, boundary_layer_top]
    maximum_loss_rate = zero(reference)

    for k in boundary_layer_top:nlayers
        loss_rate = zero(reference)
        if k > boundary_layer_top
            loss_rate += ∇²_above[k] *
                         (sigma_diffusivity[ij, k] +
                          sigma_diffusivity[ij, k - 1])
        end
        if k < nlayers
            loss_rate += ∇²_below[k] *
                         (sigma_diffusivity[ij, k] +
                          sigma_diffusivity[ij, k + 1])
        end
        maximum_loss_rate = max(maximum_loss_rate, loss_rate)
    end

    courant = maximum_timestep_seconds * maximum_loss_rate
    limit = oftype(courant, METRIC_VERTICAL_DIFFUSION_COURANT_LIMIT)
    scale = courant > limit ? limit / courant : one(courant)
    for k in boundary_layer_top:nlayers
        sigma_diffusivity[ij, k] *= scale
    end
    return scale
end

Base.@propagate_inbounds function SpeedyWeather.parameterization!(
    ij,
    vars,
    corrected::MetricBulkRichardsonDiffusion,
    model,
)
    diffusion = corrected.diffusion
    (; diffuse_momentum, diffuse_static_energy, diffuse_humidity) = diffusion
    any((diffuse_momentum, diffuse_static_energy, diffuse_humidity)) ||
        return nothing

    sigma_diffusivity, boundary_layer_top =
        SpeedyWeather.get_diffusion_coefficients!(
            ij,
            vars,
            diffusion,
            model.atmosphere,
            model.planet,
            model.orography,
            model.geopotential,
        )

    # Parameterized physics supplies a tendency to SpeedyWeather's leapfrog
    # update A_new = A_old + 2Δt * tendency.  A dissipative operator evaluated
    # on the middle state has the unstable two-level leapfrog diffusion mode.
    # Evaluating from A_old (the `_prev` fields, as other SpeedyWeather physics
    # does) instead gives a conventional one-level explicit diffusion step.
    temperature = vars.grid.temperature_prev
    humidity = vars.grid.humidity_prev
    u = vars.grid.u_prev
    v = vars.grid.v_prev
    sigma = model.geometry.σ_levels_full
    gravity = model.planet.gravity
    gas_constant = model.atmosphere.R_dry
    for k in 1:length(sigma)
        virtual_temperature = SpeedyWeather.virtual_temperature(
            temperature[ij, k],
            humidity[ij, k],
            model.atmosphere,
        )
        sigma_per_height =
            sigma[k] * gravity / (gas_constant * virtual_temperature)
        sigma_diffusivity[ij, k] *= sigma_per_height^2
    end

    # Later leapfrog calls advance across 2Δt.  This is also a safe upper bound
    # for the shorter startup steps, so the limiter does not need mutable clock
    # state inside the GPU column kernel.
    maximum_timestep_seconds = 2 * model.time_stepping.Δt_sec
    limiter_scale = _limit_metric_sigma_diffusivity!(
        ij,
        sigma_diffusivity,
        boundary_layer_top,
        diffusion,
        maximum_timestep_seconds,
    )
    vars.parameterizations.vertical_diffusion_cfl_scale[ij] = limiter_scale

    diffuse_momentum && _metric_vertical_diffusion!(
        ij,
        vars.tendencies.grid.u,
        u,
        sigma_diffusivity,
        boundary_layer_top,
        diffusion,
    )
    diffuse_momentum && _metric_vertical_diffusion!(
        ij,
        vars.tendencies.grid.v,
        v,
        sigma_diffusivity,
        boundary_layer_top,
        diffusion,
    )
    if model.atmosphere isa SpeedyWeather.AbstractWetAtmosphere &&
            diffuse_humidity
        _metric_vertical_diffusion!(
            ij,
            vars.tendencies.grid.humidity,
            humidity,
            sigma_diffusivity,
            boundary_layer_top,
            diffusion,
        )
    end

    if diffuse_static_energy
        dry_static_energy = vars.scratch.grid.a
        heat_capacity = model.atmosphere.heat_capacity
        geopotential = vars.grid.geopotential
        for k in 1:length(sigma)
            dry_static_energy[ij, k] =
                heat_capacity * temperature[ij, k] + geopotential[ij, k]
            sigma_diffusivity[ij, k] /= heat_capacity
        end
        _metric_vertical_diffusion!(
            ij,
            vars.tendencies.grid.temperature,
            dry_static_energy,
            sigma_diffusivity,
            boundary_layer_top,
            diffusion,
        )
    end
    return nothing
end
