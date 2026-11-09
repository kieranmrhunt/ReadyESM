"""v6 (6 Nov 2026, Paul's runtime-provenance requirement): the receipt and the prepared file carry a runtime
fingerprint (runtime_fingerprint): Python implementation and full version, the numpy, scipy, netCDF4 (with its netCDF-C
and HDF5 library versions) and mpmath versions, mpmath's backend, and for each of those four packages a SHA-256 over
every .py and shared-library file of its installed tree (sorted relative paths and file hashes), plus their combined
digest, the interpreter binary's SHA-256 and the environment's netCDF-C and HDF5 shared libraries' SHA-256.
verify_prepared_v6.py requires the running interpreter's fingerprint to equal the receipt's, the file's and a
pinned digest exactly. Numerics are unchanged from v5. The verdict line is UOEXETER_PREPARE_V6_DONE. The v5 text
follows.

v5 (6 Nov 2026, Paul's finding on verifier job 65711449): every number written to the output or the receipt is now
computed deterministically, independent of CPU, BLAS, LAPACK and SIMD code paths (v4 gave last-bit differences,
max 6.6e-16 relative, between a login-node and a compute-node run):
  - 1993-mean ssi: an explicitly ordered element-wise sum of the twelve months, then / 12;
  - SW weights: overlap length x value per solar bin (element-wise), summed with math.fsum (correctly rounded);
  - LW weights: Planck at 250 K integrated with mpmath (pure-Python backend, 30 digits, tanh-sinh), then rounded to
    float64 (v4: 64-point Gauss-Legendre nodes from a LAPACK eigen-solve and numpy's vectorised expm1);
  - remap: an explicitly ordered accumulation over source wavelengths, out = out + M_k * W[k, :] in source order
    (element-wise IEEE multiply and add only), then / F_b (v4: a BLAS tensor contraction);
  - covered fractions: math.fsum over each band's weights, / F_b;
  - the self-test still recomputes every weight by an independent path (scipy quad, numpy) and must pass within 1e-9,
    but only its pass flags and tolerance enter the receipt and the file; the error magnitudes, which depend on
    libm and SIMD, are printed only.
The verdict line is UOEXETER_PREPARE_V5_DONE. The v4 text follows.

v4 (6 Nov 2026, Paul's review of v3 5d8fb04e): the cross-file identity check compares every coordinate centre AND
bounds array (lat, lat_bnds, height, height_bnds, wavelength, wavelength_bnds, time, time_bnds) of ssa and asy with
ext; v3 compared lat, height and time centres only (plus wavelength bounds), so files with equal centres but different
latitude, height or time cell bounds passed. The verdict line is UOEXETER_PREPARE_V4_DONE. The v3 text follows.

v3 (6 Nov 2026, Paul's review of v2 1d0258db): (a) the coordinate BOUNDS values are validated in every source file
before extraction (validate_axis): finite, shape (n, 2), lower < upper, centres strictly increasing and each inside
its own bounds, adjacent bounds contiguous (relative 1e-12), and the axis inside its physical range (lat [-90, 90],
height [0, 100 km], wavelength (0, 1 cm), time: monthly bounds containing each centre); (b) the receipt is written to
a named RECEIPT file (not standard output) with the raw, solar and RRTMGP hashes, the script hash and the output hash,
and the prepared file carries the same input and script hashes as attributes, so verify_prepared_v3.py can bind raw
files -> receipt -> prepared output; (c) the verdict line is UOEXETER_PREPARE_V3_DONE. The v2 text follows.

v2 (6 Nov 2026, Paul's review of v1 35f4ea8a): (i) a valid cell with g < 0 or g >= 1 (or ssa outside [0, 1], or
ext < 0) fails the preparation instead of being counted, and the emitted band moments must satisfy E >= 0,
0 <= S <= E and 0 <= G < S where S > 0 (G = 0 where S = 0); (ii) before any extraction every source variable is
checked by dimension NAME (exactly time, lat, height, wavelength, in any file order, then reordered to that order),
its units (ext m-1, ssa and asy 1), the absence of fill attributes, and the coordinates (lat degrees_north, height m,
wavelength m, time days since 1850-01-01 proleptic_gregorian, bounds as (coordinate, bnds)), all equal across the
three files. The v1 text follows.

Prepare the CMIP7 stratospheric aerosol optics (UOEXETER-CMIP-2-2-1) for ReadyESM PA2 component C, December 1992 to
January 1994, on RRTMGP bands (Claude Code, 6 Nov 2026; design v4 3a0a84d1, Kieran 14:0x chose this source).

Missing data (Paul 14:4x, the only approved contract): the source marks altitudes it treats as tropospheric with NaN
(PCMDI/input4MIPs_CVs discussion #175). For every source cell (month, lat, height, wavelength) the moments
E = ext, S = ext ssa, G = ext ssa g are formed; where any of ext, ssa, g is missing all three moments are set to zero
before any remap, and those cells are counted. Any other invalid value (negative ext, ssa outside [0, 1], g outside
[-1, 1]) is refused.
Spectral remap (ReadyESM's own, conservative): each source wavelength bin [lo, hi] overlaps RRTMGP band b (limits
from rrtmgp-data v1.9, pinned hashes) over a wavelength interval; its weight W[k, b] is the integral over that overlap
of a declared spectral weight F: for SW bands the CMIP6 solar spectral irradiance (SOLARIS-HEPPA-3-2 ssi, 1993 mean,
piecewise constant on its own bins, integrated exactly through its cumulative integral), for LW bands the Planck
function at 250 K (64-point Gauss-Legendre per overlap). Band moment = sum_k W[k, b] M_k / F_b, where F_b is the
integral of F over the WHOLE band: band parts outside the source range (125 nm to 112.5 um; only LW band 1, 10-250
cm-1, is partly outside) contribute zero, as missing cells do. ssa = S / E and g = G / S at use, only where the
denominators are nonzero. The weight table and each band's covered fraction are written to the output.
Self-tests (independent code paths): (1) the solar overlap integrals recomputed by a direct loop over the solar bins
(overlap length x value) and the Planck integrals by scipy.integrate.quad must equal W (and the whole-band weights) within 1e-9 relative;
(2) a constant source field maps to constant x covered fraction in every band.
usage: python3 -I prepare_uoexeter_v6.py UOEXETER_DIR SOLAR_FILE RRTMGP_DATA_DIR OUTPUT.nc RECEIPT.json
"""
import hashlib
import json
import math
import os
import platform
import sys

import netCDF4
import numpy as np
import mpmath
from scipy import integrate

SOURCES = {
    "ext": ("ext_input4MIPs_aerosolProperties_CMIP_UOEXETER-CMIP-2-2-1_gnz_175001-202312.nc",
            "2a3980ff420487f94816f53f3fa2b6518bd28104d30ad66d55dd5f24baabc422"),
    "ssa": ("ssa_input4MIPs_aerosolProperties_CMIP_UOEXETER-CMIP-2-2-1_gnz_175001-202312.nc",
            "c7406b8651a629eaec6cf9d84614e6b838ba5ff0d0b948dc60bdc7000b6f5abe"),
    "asy": ("asy_input4MIPs_aerosolProperties_CMIP_UOEXETER-CMIP-2-2-1_gnz_175001-202312.nc",
            "1f3f3f984b5c788cabe8df682f7ff31205688c4733575ee198339976e9154f60"),
}
SOLAR_SHA256 = "4b844a3324b74be16faf7272b64ed7853ea2c606c4db841a318366ad223ebcac"
RRTMGP_SHA256 = {"rrtmgp-gas-sw-g224.nc": "584f1dd41ea9fc07d4ee3754eb1dafbd46ad3161cd6fd20fa06b6922b6f0702e",
                 "rrtmgp-gas-lw-g256.nc": "4048360199d1917ed8f2ccaae2ec097d0f990da3bbad9830337b739b4fa01be7"}
MONTHS = [(1992, 12)] + [(1993, m) for m in range(1, 13)] + [(1994, 1)]
PLANCK_T = 250.0
H, C, KB = 6.62607015e-34, 2.99792458e8, 1.380649e-23


def _tree_digest(package):
    """SHA-256 over every .py and shared-library file of an installed package tree: sorted relative paths and hashes."""
    root = os.path.dirname(os.path.abspath(package.__file__))
    lines = []
    for d, dirs, files in os.walk(root):
        dirs.sort()
        for f in sorted(files):
            if f.endswith((".py", ".so")) or ".so." in f:
                p = os.path.join(d, f)
                lines.append(os.path.relpath(p, root) + " " + sha256(p))
    return hashlib.sha256("\n".join(sorted(lines)).encode()).hexdigest(), len(lines)


def runtime_fingerprint():
    """The numerical runtime this preparation depends on (v6, Paul): interpreter and library identities."""
    import scipy
    fp = {"python_implementation": platform.python_implementation(), "python_version": sys.version,
          "numpy": np.__version__, "scipy": scipy.__version__, "netCDF4": netCDF4.__version__,
          "netcdf_c": netCDF4.__netcdf4libversion__, "hdf5": netCDF4.__hdf5libversion__,
          "mpmath": mpmath.__version__, "mpmath_backend": mpmath.libmp.BACKEND}
    for name, mod in (("numpy", np), ("scipy", scipy), ("netCDF4", netCDF4), ("mpmath", mpmath)):
        fp[f"{name}_tree_sha256"], fp[f"{name}_tree_files"] = _tree_digest(mod)
    # The interpreter binary and the netCDF-C and HDF5 shared libraries of this environment (not inside the trees).
    fp["python_executable_sha256"] = sha256(os.path.realpath(sys.executable))
    lib = os.path.join(sys.prefix, "lib")
    shared = sorted(f for f in os.listdir(lib) if f.startswith(("libnetcdf.so", "libhdf5.so", "libhdf5_hl.so"))
                    and not os.path.islink(os.path.join(lib, f)))
    if not shared:
        raise ValueError("netCDF-C / HDF5 shared libraries not found in " + lib)
    fp["shared_libraries_sha256"] = {f: sha256(os.path.join(lib, f)) for f in shared}
    fp["digest"] = hashlib.sha256(json.dumps(fp, sort_keys=True).encode()).hexdigest()
    return fp


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


def validate_axis(label, name, centre, bnds, lo_limit, hi_limit, contiguous=True):
    """Refuse malformed coordinate bounds (v3, Paul): every check raises ValueError naming the file and axis."""
    centre, bnds = np.asarray(centre, dtype=float), np.asarray(bnds, dtype=float)
    def fail(why):
        raise ValueError(f"{label}: {name} bounds {why}")
    if bnds.ndim != 2 or bnds.shape != (centre.size, 2):
        fail(f"have shape {bnds.shape}, expected ({centre.size}, 2)")
    if not (np.all(np.isfinite(bnds)) and np.all(np.isfinite(centre))):
        fail("or centres are not finite")
    lo, hi = bnds[:, 0], bnds[:, 1]
    if not np.all(lo < hi):
        fail("have a lower bound not below its upper bound")
    if not np.all(np.diff(centre) > 0):
        fail("belong to centres that are not strictly increasing")
    if not np.all((lo <= centre) & (centre <= hi)):
        fail("do not contain their centres")
    if contiguous and not np.allclose(lo[1:], hi[:-1], rtol=1e-12, atol=0):
        fail("are not contiguous")
    if not (lo[0] >= lo_limit and hi[-1] <= hi_limit):
        fail(f"span {lo[0]}..{hi[-1]}, outside {lo_limit}..{hi_limit}")
    return {"n": int(centre.size), "lower": float(lo[0]), "upper": float(hi[-1])}


def planck(lam_m, t=PLANCK_T):
    return 2 * H * C**2 / lam_m**5 / np.expm1(H * C / (lam_m * KB * t))


def solar_bins(path):
    """1993-mean ssi (W m-2 nm-1) on its own bins, bounds in m, ascending."""
    with netCDF4.Dataset(path) as ds:
        ds.set_auto_mask(False)
        rows = np.where(ds["calyear"][:] == 1993)[0]
        assert len(rows) == 12, "solar file lacks the twelve 1993 months"
        raw = ds["ssi"][rows, :].astype(float)
        ssi = np.zeros(raw.shape[1])
        for r in range(raw.shape[0]):                                     # explicitly ordered, element-wise
            ssi = ssi + raw[r]
        ssi = ssi / 12.0
        bnds = ds["wlen_bnds"][:].astype(float) * 1e-9
    order = np.argsort(bnds[:, 0])
    lo, hi, val = bnds[order, 0], bnds[order, 1], ssi[order]
    assert np.allclose(lo[1:], hi[:-1], rtol=1e-12) and np.all(np.isfinite(val)) and np.all(val >= 0)
    return lo, hi, val


def solar_integral(bins):
    """Integral of the piecewise-constant spectrum from a to b (zero outside): overlap length x value per solar bin,
    element-wise, summed with math.fsum (correctly rounded, so independent of summation order and SIMD)."""
    lo, hi, val = bins
    def f(a, b):
        ov = np.clip(np.minimum(hi, b) - np.maximum(lo, a), 0.0, None)
        return math.fsum((ov * val).tolist())
    return f


_MP = mpmath.mp
_MP.dps = 30
_H, _C, _KB, _T = (_MP.mpf(x) for x in ("6.62607015e-34", "2.99792458e8", "1.380649e-23", "250"))


def planck_integral(a, b):
    """Planck (250 K) integrated over [a, b] m with mpmath at 30 digits (pure Python, deterministic), as float64."""
    f = lambda lam: 2 * _H * _C**2 / lam**5 / _MP.expm1(_H * _C / (lam * _KB * _T))
    return float(_MP.quad(f, [_MP.mpf(a), _MP.mpf(b)]))


def band_limits(rrtmgp, name):
    path = os.path.join(rrtmgp, name)
    assert sha256(path) == RRTMGP_SHA256[name], name
    with netCDF4.Dataset(path) as ds:
        wn = ds["bnd_limits_wavenumber"][:].astype(float)               # cm-1, (band, 2)
    return [(0.01 / b[1], 0.01 / b[0]) for b in wn]                     # wavelength m (lo, hi)


def weights(src_lo, src_hi, bands, integral):
    """W[k, b] = integral of F over the overlap of source bin k and band b; Fb[b] = integral of F over band b."""
    W = np.zeros((len(src_lo), len(bands)))
    Fb = np.array([integral(blo, bhi) for blo, bhi in bands])
    for b, (blo, bhi) in enumerate(bands):
        for k in range(len(src_lo)):
            lo, hi = max(src_lo[k], blo), min(src_hi[k], bhi)
            if hi > lo:
                W[k, b] = integral(lo, hi)
    assert np.all(Fb > 0)
    return W, Fb


def remap(moments, W, Fb):
    """Band moments = sum_k W[k, b] M_k / Fb[b]; last axis of `moments` is the source wavelength. Accumulated in source
    order with element-wise IEEE multiply and add only (no BLAS, no reduction reordering)."""
    moments = np.asarray(moments, dtype=np.float64)
    out = np.zeros(moments.shape[:-1] + (W.shape[1],))
    for k in range(W.shape[0]):
        out = out + moments[..., k, None] * W[k, :]
    return out / Fb


def covered_fraction(W, Fb):
    """Each band's covered fraction, math.fsum of its weights / F_b (deterministic)."""
    return [math.fsum(W[:, b].tolist()) / Fb[b] for b in range(W.shape[1])]


def self_test(src_lo, src_hi, bands, W, Fb, kind, sbins=None):
    ref = np.zeros_like(W)
    for b, (blo, bhi) in enumerate(bands):
        for k in range(len(src_lo)):
            lo, hi = max(src_lo[k], blo), min(src_hi[k], bhi)
            if hi <= lo:
                continue
            if kind == "sw":
                slo, shi, val = sbins
                ov = np.clip(np.minimum(shi, hi) - np.maximum(slo, lo), 0.0, None)
                ref[k, b] = float(np.sum(ov * val))
            else:
                ref[k, b] = integrate.quad(planck, lo, hi, epsrel=1e-12, limit=200)[0]
    if kind == "sw":
        slo, shi, val = sbins
        refF = np.array([float(np.sum(np.clip(np.minimum(shi, bhi) - np.maximum(slo, blo), 0.0, None) * val))
                         for blo, bhi in bands])
    else:
        refF = np.array([integrate.quad(planck, blo, bhi, epsrel=1e-12, limit=400)[0] for blo, bhi in bands])
    errF = float(np.max(np.abs(Fb - refF) / refF))
    nz = ref > 0
    err = max(float(np.max(np.abs(W[nz] - ref[nz]) / ref[nz])), errF)
    same_support = bool(np.array_equal(nz, W > 0))
    frac = np.array(covered_fraction(W, Fb))
    const = remap(np.full(len(src_lo), 3.5), W, Fb)
    ok_const = bool(np.allclose(const, 3.5 * frac, rtol=1e-12))
    return (err <= 1e-9 and same_support and ok_const), err, [float(f) for f in frac]


def main(src_dir, solar_path, rrtmgp, output, receipt_path):
    if os.path.exists(output) or os.path.exists(receipt_path):
        sys.exit("refusing to overwrite " + output + " or " + receipt_path)
    receipt = {"inputs": {}, "script_sha256": sha256(os.path.abspath(__file__)), "runtime": runtime_fingerprint()}
    for v, (name, pinned) in SOURCES.items():
        h = sha256(os.path.join(src_dir, name))
        assert h == pinned, (name, h)
        receipt["inputs"][name] = h
        receipt.setdefault("input_bytes", {})[name] = os.path.getsize(os.path.join(src_dir, name))
    for name in RRTMGP_SHA256:
        receipt["inputs"][name] = sha256(os.path.join(rrtmgp, name))
    hs = sha256(solar_path)
    assert hs == SOLAR_SHA256, hs
    receipt["inputs"][os.path.basename(solar_path)] = hs
    data, coords = {}, {}
    order = ("time", "lat", "height", "wavelength")
    units = {"ext": "m-1", "ssa": "1", "asy": "1"}
    coord_units = {"lat": "degrees_north", "height": "m", "wavelength": "m"}
    for v, (name, _) in SOURCES.items():
        with netCDF4.Dataset(os.path.join(src_dir, name)) as ds:
            ds.set_auto_mask(False)
            var = ds[v]
            dims = var.dimensions
            if sorted(dims) != sorted(order) or len(dims) != 4:
                raise ValueError(f"{name}: {v} has dimensions {dims}, expected {order} in some order")
            if getattr(var, "units", None) != units[v]:
                raise ValueError(f"{name}: {v} units {getattr(var, 'units', None)!r}, expected {units[v]!r}")
            for a in ("_FillValue", "missing_value"):
                if a in var.ncattrs():
                    raise ValueError(f"{name}: {v} declares {a}")
            for c, u in coord_units.items():
                if ds[c].dimensions != (c,) or getattr(ds[c], "units", None) != u:
                    raise ValueError(f"{name}: coordinate {c} has dimensions {ds[c].dimensions} units "
                                     f"{getattr(ds[c], 'units', None)!r}, expected ({c},) {u!r}")
                b = ds[c + "_bnds"] if c != "wavelength" else ds["wavelength_bnds"]
                if b.dimensions[0] != c or ds.dimensions[b.dimensions[1]].size != 2:
                    raise ValueError(f"{name}: {c} bounds have dimensions {b.dimensions}")
            if ds["time"].units != "days since 1850-01-01" or ds["time"].calendar != "proleptic_gregorian":
                raise ValueError(f"{name}: time is {ds['time'].units} {ds['time'].calendar}")
            if ds["time_bnds"].dimensions[0] != "time":
                raise ValueError(f"{name}: time bounds have dimensions {ds['time_bnds'].dimensions}")
            axes = receipt.setdefault("axis_bounds", {})
            for c, b, lo_l, hi_l, contig in (("lat", "lat_bnds", -90.0, 90.0, True),
                                             ("height", "height_bnds", 0.0, 1.0e5, True),
                                             ("wavelength", "wavelength_bnds", 0.0, 1.0e-2, True),
                                             ("time", "time_bnds", -1.0e7, 1.0e7, False)):
                axes[f"{v}.{c}"] = validate_axis(name, c, ds[c][:], ds[b][:], lo_l, hi_l, contig)
            if receipt["axis_bounds"][f"{v}.wavelength"]["lower"] <= 0:
                raise ValueError(f"{name}: wavelength bounds are not positive")
            coords[v] = [ds[c][:].astype(float) for c in ("lat", "lat_bnds", "height", "height_bnds", "wavelength",
                                                          "wavelength_bnds", "time", "time_bnds")]
            t = netCDF4.num2date(ds["time"][:], ds["time"].units, ds["time"].calendar)
            idx = [next(i for i, d in enumerate(t) if (d.year, d.month) == ym) for ym in MONTHS]
            sel = tuple(idx if d == "time" else slice(None) for d in dims)
            raw = np.asarray(var[sel], dtype=np.float64)
            data[v] = np.transpose(raw, [dims.index(d) for d in order])       # (time, lat, height, wavelength)
            if v == "ext":
                lat = ds["lat"][:].astype(float); lat_b = ds["lat_bnds"][:].astype(float)
                hgt = ds["height"][:].astype(float); hgt_b = ds["height_bnds"][:].astype(float)
                wb = ds["wavelength_bnds"][:].astype(float)
                licence = ds.getncattr("license"); src_id = ds.getncattr("source_id")
                tb = ds["time_bnds"][idx, :]
                tunits, tcal = ds["time"].units, ds["time"].calendar
    receipt["source_dimensions_validated"] = {"order": list(order), "units": units, "coordinates": coord_units}
    for v in ("ssa", "asy"):
        names = ("lat", "lat_bnds", "height", "height_bnds", "wavelength", "wavelength_bnds", "time", "time_bnds")
        differ = [n for n, a, b in zip(names, coords["ext"], coords[v]) if not np.array_equal(a, b)]
        if differ:
            raise ValueError(f"{v} coordinates differ from ext: {differ}")
    ext, ssa, asy = data["ext"], data["ssa"], data["asy"]
    missing = ~(np.isfinite(ext) & np.isfinite(ssa) & np.isfinite(asy))
    valid = ~missing
    # ReadyESM admits 0 <= ssa <= 1 and 0 <= g < 1 only (Paul): any valid cell outside fails the preparation.
    bad = valid & ((ext < 0) | (ssa < 0) | (ssa > 1) | (asy < 0) | (asy >= 1))
    if bad.any():
        raise ValueError(f"{int(bad.sum())} valid source cells have optics outside ext >= 0, 0 <= ssa <= 1, "
                         f"0 <= g < 1 (negative g: {int((valid & (asy < 0)).sum())}, g >= 1: "
                         f"{int((valid & (asy >= 1)).sum())})")
    E = np.where(valid, ext, 0.0)
    S = np.where(valid, ext * ssa, 0.0)
    G = np.where(valid, ext * ssa * asy, 0.0)
    receipt["zeroed_missing_cells"] = int(missing.sum())
    # Diagnostic only: the provider's NaN cells are tropospheric, so in each column they should lie below every valid
    # cell. Columns where a missing cell sits above a valid one are counted, not refused.
    col_missing = missing.all(axis=-1)                                   # (time, lat, height)
    above_valid = np.zeros(col_missing.shape[:2], dtype=bool)
    for k in range(1, col_missing.shape[2]):
        above_valid |= col_missing[:, :, k] & ~col_missing[:, :, :k].all(axis=2)
    receipt["columns_with_missing_above_valid"] = int(above_valid.sum())
    receipt["missing_mixed_over_wavelength"] = int((missing.any(axis=-1) & ~missing.all(axis=-1)).sum())
    receipt["source_cells"] = int(missing.size)
    order = np.argsort(wb[:, 0])
    src_lo, src_hi = wb[order, 0], wb[order, 1]
    E, S, G = E[..., order], S[..., order], G[..., order]
    sbins = solar_bins(solar_path)
    sw = band_limits(rrtmgp, "rrtmgp-gas-sw-g224.nc")
    lw = band_limits(rrtmgp, "rrtmgp-gas-lw-g256.nc")
    W_sw, F_sw = weights(src_lo, src_hi, sw, solar_integral(sbins))
    W_lw, F_lw = weights(src_lo, src_hi, lw, planck_integral)
    out = {}
    for spec, W, Fb in (("sw", W_sw, F_sw), ("lw", W_lw, F_lw)):
        for nm, M in (("E", E), ("S", S), ("G", G)):
            out[f"{nm}_{spec}"] = remap(M, W, Fb)
    ok_sw, err_sw, frac_sw = self_test(src_lo, src_hi, sw, W_sw, F_sw, "sw", sbins)
    ok_lw, err_lw, frac_lw = self_test(src_lo, src_hi, lw, W_lw, F_lw, "lw")
    receipt["sw_covered_fraction"], receipt["lw_covered_fraction"] = frac_sw, frac_lw
    # Only deterministic pass flags enter the receipt and the file; the error magnitudes depend on libm/SIMD paths.
    receipt["self_test"] = {"sw_pass": ok_sw, "lw_pass": ok_lw, "tolerance_rel": 1e-9}
    print(f"SELF_TEST_DIAGNOSTIC sw_max_rel_err={err_sw:.3e} lw_max_rel_err={err_lw:.3e}")
    if not (ok_sw and ok_lw):
        raise ValueError(f"spectral remap self-test failed: {receipt['self_test']}")
    for k, v in out.items():
        if not np.all(np.isfinite(v)) or (not k.startswith("G") and np.any(v < 0)):
            raise ValueError(f"remapped {k} has nonfinite or negative values")
    for spec in ("sw", "lw"):
        E_, S_, G_ = out[f"E_{spec}"], out[f"S_{spec}"], out[f"G_{spec}"]
        if np.any(S_ > E_ * (1 + 1e-12)) or np.any(G_ < 0) or np.any((S_ > 0) & (G_ >= S_)) or \
                np.any((S_ == 0) & (G_ != 0)):
            raise ValueError(f"remapped {spec} moments violate S <= E, 0 <= G < S (G = 0 where S = 0)")
    with netCDF4.Dataset(output, "w", format="NETCDF4") as ds:
        ds.schema = "readiesm_strat_aerosol_bands_v1"
        ds.source = (f"{src_id} (CMIP7 input4MIPs, esgf-node.ornl.gov, v20250521); provider code "
                     "thomasaubry/CMIP7_stratforcing_v2.2.1 head f2098b7c (provenance only)")
        ds.licence = licence
        ds.input_sha256 = json.dumps(receipt["inputs"])
        ds.missing_rule = ("source NaN cells (tropospheric per the provider, input4MIPs_CVs #175): E = ext, S = ext ssa, "
                           "G = ext ssa g all set to zero before remap; ssa = S/E, g = G/S at use where nonzero")
        ds.spectral_rule = ("conservative overlap of source wavelength bins with RRTMGP bands (rrtmgp-data v1.9), "
                            "weighted by the 1993-mean CMIP6 SOLARIS-HEPPA-3-2 ssi (SW) or Planck 250 K (LW); "
                            "moment-preserving (E, S, G); normalised by the whole-band weight, so band parts outside "
                            "the source range contribute zero")
        ds.label = "UOEXETER-CMIP-2-2-1 (CMIP7) stratospheric optics, observation-based, for PA2 component C"
        ds.zeroed_missing_cells = receipt["zeroed_missing_cells"]
        ds.self_test = json.dumps(receipt["self_test"])
        ds.preparation_script_sha256 = receipt["script_sha256"]
        ds.runtime_fingerprint = json.dumps(receipt["runtime"], sort_keys=True)
        for nm, n in (("time", len(MONTHS)), ("lat", len(lat)), ("height", len(hgt)), ("sw_band", len(sw)),
                      ("lw_band", len(lw)), ("source_wavelength", len(src_lo)), ("bnds", 2)):
            ds.createDimension(nm, n)
        v = ds.createVariable("time_bnds", "f8", ("time", "bnds"), fill_value=False); v[:] = tb
        v.units = tunits; v.calendar = tcal
        ds.createVariable("year", "i4", ("time",), fill_value=False)[:] = [y for y, _ in MONTHS]
        ds.createVariable("month", "i4", ("time",), fill_value=False)[:] = [m for _, m in MONTHS]
        ds.createVariable("lat", "f8", ("lat",), fill_value=False)[:] = lat
        ds.createVariable("lat_bnds", "f8", ("lat", "bnds"), fill_value=False)[:] = lat_b
        v = ds.createVariable("height", "f8", ("height",), fill_value=False); v[:] = hgt; v.units = "m above sea level"
        ds.createVariable("height_bnds", "f8", ("height", "bnds"), fill_value=False)[:] = hgt_b
        for k, arr in out.items():
            band = "sw_band" if k.endswith("sw") else "lw_band"
            v = ds.createVariable(k, "f8", ("time", "lat", "height", band), fill_value=False)
            v[:] = arr
            v.units = "m-1"
        ds.createVariable("weights_sw", "f8", ("source_wavelength", "sw_band"), fill_value=False)[:] = W_sw
        ds.createVariable("weights_lw", "f8", ("source_wavelength", "lw_band"), fill_value=False)[:] = W_lw
        ds.createVariable("band_weight_sw", "f8", ("sw_band",), fill_value=False)[:] = F_sw
        ds.createVariable("band_weight_lw", "f8", ("lw_band",), fill_value=False)[:] = F_lw
        ds.createVariable("covered_fraction_sw", "f8", ("sw_band",), fill_value=False)[:] = frac_sw
        ds.createVariable("covered_fraction_lw", "f8", ("lw_band",), fill_value=False)[:] = frac_lw
        ds.createVariable("source_wavelength_lo_m", "f8", ("source_wavelength",), fill_value=False)[:] = src_lo
        ds.createVariable("source_wavelength_hi_m", "f8", ("source_wavelength",), fill_value=False)[:] = src_hi
    receipt["output_sha256"] = sha256(output)
    receipt["verdict"] = "UOEXETER_PREPARE_V6_DONE"
    with open(receipt_path, "x") as f:
        json.dump(receipt, f, indent=1, sort_keys=True)
    print("UOEXETER_PREPARE_V6_DONE output_sha256=" + receipt["output_sha256"] + " runtime=" + receipt["runtime"]["digest"])


if __name__ == "__main__":
    if len(sys.argv) != 6:
        sys.exit(__doc__)
    main(*sys.argv[1:])
