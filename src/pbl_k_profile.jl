#=
Boundary-layer K-profile mixing: column functions (Stage 1 prototype, not yet connected to the model).

Production has no vertical diffusion, so Betts-Miller convection is the only process that removes heat and moisture
from the lower troposphere. This file holds the column physics of an opt-in replacement, after Holtslag and Boville
(1993), first order and without a counter-gradient term:

- the boundary-layer depth h is the first height where the bulk Richardson number reaches the critical value;
- inside the boundary layer K(z) = kappa w_s z (1 - z/h)^2 with the velocity scale w_s from u* and, for an unstable
  surface buoyancy flux, the convective velocity w*; above h, K = 0;
- dry static energy and specific humidity are mixed by a backward-Euler tridiagonal solve in sigma flux form with no
  flux at the surface and at the top, which is unconditionally stable and conserves the sigma-weighted column sum.

Index convention as in SpeedyWeather: k = 1 is the top layer, k = n the lowest. Half level k sits above full level k,
so half level k + 1 is the interface between layers k and k + 1, and half level n + 1 is the surface.
=#

const PBL_VON_KARMAN = 0.4
const PBL_CRITICAL_RICHARDSON = 0.3
const PBL_THERMAL_EXCESS_FACTOR = 8.5
const PBL_CONVECTIVE_VELOCITY_FACTOR = 0.6
const PBL_STABLE_PROFILE_FACTOR = 5.0
const PBL_SHEAR_FLOOR_FACTOR = 100.0
# Calm limit: with zero resolved wind and zero u* the bulk Richardson number is 0/0. A 0.1 m/s wind floor in its
# denominator keeps it finite (stable calm columns get a large Ri and a shallow layer, neutral ones Ri = 0). It
# changes Ri by about 0.04 % at a typical 5 m/s wind.
const PBL_CALM_WIND_FLOOR = 0.1

"""Obukhov length (m) for friction velocity `ustar` and surface kinematic virtual heat flux `buoyancy_flux` (K m/s).
Infinite for zero flux, positive for a stable (downward) flux, negative for an unstable one."""
function pbl_obukhov_length(ustar, buoyancy_flux, thetav_reference, gravity)
    buoyancy_flux == 0 && return oftype(float(ustar), Inf)
    return -ustar^3 * thetav_reference / (PBL_VON_KARMAN * gravity * buoyancy_flux)
end

"""Convective velocity scale w* (m/s); zero unless the surface buoyancy flux is upward."""
pbl_convective_velocity(h, buoyancy_flux, thetav_reference, gravity) =
    buoyancy_flux > 0 ? cbrt(gravity / thetav_reference * buoyancy_flux * h) : zero(float(h))

"""Turbulent velocity scale w_s (m/s) at height `z` in a boundary layer of depth `h`."""
function pbl_velocity_scale(z, h, ustar, buoyancy_flux, thetav_reference, gravity)
    # Calm and stable or neutral: no turbulence velocity (avoids 0/0 in the Obukhov length).
    ustar == 0 && buoyancy_flux <= 0 && return zero(float(ustar))
    if buoyancy_flux > 0
        wstar = pbl_convective_velocity(h, buoyancy_flux, thetav_reference, gravity)
        return cbrt(ustar^3 + PBL_CONVECTIVE_VELOCITY_FACTOR * wstar^3)
    end
    obukhov = pbl_obukhov_length(ustar, buoyancy_flux, thetav_reference, gravity)
    return ustar / (1 + PBL_STABLE_PROFILE_FACTOR * z / obukhov)
end

"""Eddy diffusivity K (m2/s) at height `z` above the surface; zero at and above `h`."""
function pbl_diffusivity(z, h, ustar, buoyancy_flux, thetav_reference, gravity)
    (0 < z < h) || return zero(float(z))
    w = pbl_velocity_scale(z, h, ustar, buoyancy_flux, thetav_reference, gravity)
    return PBL_VON_KARMAN * w * z * (1 - z / h)^2
end

"""Height (m) where the bulk Richardson number first reaches `critical` scanning up from the lowest full level,
linearly interpolated, for a surface reference virtual potential temperature `thetav_surface`."""
function _pbl_richardson_crossing(z_full, thetav_full, wind2_full, ustar, thetav_surface, thetav_reference, gravity,
                                  critical)
    n = length(z_full)
    shear_floor = PBL_SHEAR_FLOOR_FACTOR * ustar^2 + PBL_CALM_WIND_FLOOR^2
    richardson(k) = gravity * z_full[k] * (thetav_full[k] - thetav_surface) /
                    (thetav_reference * (wind2_full[k] + shear_floor))
    previous = richardson(n)
    previous >= critical && return z_full[n]
    for k in (n - 1):-1:1
        current = richardson(k)
        if current >= critical
            fraction = (critical - previous) / (current - previous)
            return z_full[k + 1] + fraction * (z_full[k] - z_full[k + 1])
        end
        previous = current
    end
    return z_full[1]
end

"""
Boundary-layer depth (m) from the bulk Richardson criterion, clamped to [`z_min`, `z_max`]. Under an upward surface
buoyancy flux the surface parcel carries the HB93 thermal excess 8.5 B_s / w_m, with w_m from a first, excess-free
depth estimate (two passes).
"""
function pbl_height(z_full, thetav_full, wind2_full, ustar, buoyancy_flux, thetav_reference, gravity, z_min, z_max;
                    critical = PBL_CRITICAL_RICHARDSON)
    thetav_surface = thetav_full[end]
    h = _pbl_richardson_crossing(z_full, thetav_full, wind2_full, ustar, thetav_surface, thetav_reference, gravity,
                                 critical)
    if buoyancy_flux > 0
        h = clamp(h, z_min, z_max)
        wstar = pbl_convective_velocity(h, buoyancy_flux, thetav_reference, gravity)
        w_m = cbrt(ustar^3 + PBL_CONVECTIVE_VELOCITY_FACTOR * wstar^3)
        excess = w_m > 0 ? PBL_THERMAL_EXCESS_FACTOR * buoyancy_flux / w_m : zero(w_m)
        h = _pbl_richardson_crossing(z_full, thetav_full, wind2_full, ustar, thetav_surface + excess,
                                     thetav_reference, gravity, critical)
    end
    return clamp(h, z_min, z_max)
end

"""
Interface coefficients `a[k]` (1/s) for the interface between layers k and k + 1, k = 1:n-1:
a = (sigma_half g / (R_d Tv_half))^2 K(z_half) / (sigma_full[k+1] - sigma_full[k]), the hydrostatic conversion of the
height-coordinate flux rho K dx/dz to sigma. K is evaluated at the face height and is zero for faces at or above `h`.
Finite-volume consequence: a cell whose lower face lies below h exchanges through that face even when its centre is
above h (a partially intersected top cell). Cells whose upper and lower faces are both at or above h are unchanged.
"""
function pbl_interface_coefficients!(a, z_half, sigma_full, sigma_half, tv_half, h, ustar, buoyancy_flux,
                                     thetav_reference, gravity, gas_constant)
    n = length(sigma_full)
    for k in 1:(n - 1)
        z = z_half[k + 1]
        K = pbl_diffusivity(z, h, ustar, buoyancy_flux, thetav_reference, gravity)
        metric = (sigma_half[k + 1] * gravity / (gas_constant * tv_half[k + 1]))^2
        a[k] = metric * K / (sigma_full[k + 1] - sigma_full[k])
    end
    return a
end

"""
Backward-Euler step of length `dt` for dx/dt = (F_{k+1/2} - F_{k-1/2}) / dsigma_k, F_{k+1/2} = a[k] (x[k+1] - x[k]),
no flux at the top and the surface. Writes the new column into `x_new` (Thomas algorithm, `c` and `d` scratch of length
n). The matrix is a diagonally dominant M-matrix: the solve is stable, creates no new extrema, and conserves
sum dsigma x to roundoff.
"""
function pbl_implicit_diffuse!(x_new, x, a, dsigma, dt, c, d)
    n = length(x)
    # Row k: -lower x[k-1] + (1 + lower + upper) x[k] - upper x[k+1] = x_old[k].
    upper = dt * a[1] / dsigma[1]
    diagonal = 1 + upper
    c[1] = -upper / diagonal
    d[1] = x[1] / diagonal
    for k in 2:n
        lower = dt * a[k - 1] / dsigma[k]
        upper = k < n ? dt * a[k] / dsigma[k] : zero(lower)
        denominator = (1 + lower + upper) + lower * c[k - 1]
        c[k] = -upper / denominator
        d[k] = (x[k] + lower * d[k - 1]) / denominator
    end
    x_new[n] = d[n]
    for k in (n - 1):-1:1
        x_new[k] = d[k] - c[k] * x_new[k + 1]
    end
    return x_new
end

"""
Flux-form tendency (x units per second) of the backward-Euler step: solve for the new column, then form
(F_{k+1/2} - F_{k-1/2}) / dsigma_k from the new state's face fluxes. In exact arithmetic this equals (x_new - x) / dt.
In floating point it avoids subtracting two nearly equal states (about 3e5 J/kg for dry static energy), and the face
fluxes cancel pairwise in the sigma-weighted column sum.
"""
function pbl_implicit_tendency!(tendency, x, a, dsigma, dt, c, d, x_new)
    n = length(x)
    pbl_implicit_diffuse!(x_new, x, a, dsigma, dt, c, d)
    flux_above = zero(eltype(tendency))
    for k in 1:n
        flux_below = k < n ? a[k] * (x_new[k + 1] - x_new[k]) : zero(flux_above)
        tendency[k] = (flux_below - flux_above) / dsigma[k]
        flux_above = flux_below
    end
    return tendency
end
