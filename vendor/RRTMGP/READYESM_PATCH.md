# ReadyESM patch

Based on RRTMGP 0.22.2, git tree `3d9478c07e05aba9f1b97a5b94f3655f4d2e006c`, which is the tree the ReadyESM Manifest
pinned. The whole package is kept: source, extensions, tests, docs, `LICENSE` (Apache-2.0) and `LICENSE_original`.
The parent commit of this patch holds the unchanged copy.

The patch adds one optional input: prescribed aerosol optics per band. The prescribed aerosol scheme (PA2) needs it,
because it applies the extinction, single-scattering albedo and asymmetry of the CMIP aerosol data sets directly.

The hook is additive: the prescribed optics are added to whatever the MERRA lookup contributes. ReadyESM therefore
requires the MERRA aerosol mass to be zero whenever band optics are present, so there is no double counting. This is
enforced by `_enable_sulfate_aerosol(...; band_optics = true)`, which refuses a nonzero sulfate target, and by the
admission check `_validate_band_optics` (ReadyESM src/rrtmgp_radiation.jl). The check also enforces shape, device,
finite values, tau >= 0, 0 <= ssa <= 1 and 0 <= g < 1 before any solve. The kernels skip layers with tau_b = 0.
A NaN or negative tau would also be skipped silently by the kernels, which is why the check must run first.

- `AtmosphericStates.AerosolState` gains a sixth field, `band_optics`, and a fourth type parameter.
  - `nothing` keeps the unpatched behaviour.
  - Otherwise `band_optics` is a NamedTuple `(tau_sw, ssa_sw, g_sw, tau_lw, ssa_lw, g_lw)`. Each array has the shape
    `(n_bands, nlay, ncol)`.
  - A five-argument outer constructor keeps every existing caller unchanged.
- The prescribed optics are added after the MERRA lookup optics at the three places where that lookup is applied
  (`Optics.compute_optical_props!`). This is in the same way and with the same operations as the lookup:
  - **LW, no scattering:** `tau += tau_b (1 - ssa_b)`.
  - **LW, two-stream:** `increment_2stream`.
  - **SW, two-stream:** the 550 nm band's extinction and scattering go into the AOD diagnostics first. Then the optics
    are delta-scaled (`delta_scale`) and combined with `increment_2stream`.
- A layer with `tau_b = 0` is skipped. So zero prescribed optics leave the optical properties bit for bit unchanged.
  With `band_optics === nothing`, dispatch selects no-op methods.
- The prescribed optics are used only where RRTMGP already applies aerosol optics: solvers with the aerosol lookup
  loaded. The SW `OneScalar` path has no aerosol in the unpatched package and has none here.

Gates, in ReadyESM `test/rrtmgp_band_optics_hook.jl`:
- analytic one-layer values against hand-written formulas;
- bitwise flux parity of zero prescribed optics against no prescribed optics;
- the 550 nm AOD diagnostic equal to the prescribed band optical depth;
- the admission check refusing each corrupted field and a nonzero MERRA mass;
- CPU/GPU agreement, in a separate GPU job.
