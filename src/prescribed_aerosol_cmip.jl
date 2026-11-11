#=
Prescribed CMIP aerosol optics (PA2, design v3 0e5660da), opt-in with forcing.aerosol_scheme = :cmip_prescribed.

The aerosol is given to RRTMGP as optics per band (vendored RRTMGP band-optics hook, vendor/RRTMGP/READYESM_PATCH.md),
not as MERRA masses; the MERRA mass stays zero, checked at every write (`_set_band_optics!`). Components:
  A  MACv2 pre-industrial background (Kinne 2019; prepared file readiesm_macv2_background_v1): SW fine (1850) and
     coarse modes, LW total, monthly climatology on 1 degree, with 20-layer AOD fractions in height above sea level.
  B  MACv2-SP-1-0 anthropogenic plumes for the model year (macv2_simple_plumes.jl), SW only, zero above 3000 nm.
  C  stratospheric optics (CMIP7 UOEXETER-CMIP-2-2-1, prepared by prepare_uoexeter_v1.py into band moments of
     extinction per metre on the source latitudes and 500 m height layers, schema readiesm_strat_aerosol_bands_v1),
     optional (forcing.aerosol_cmip_stratospheric_path); without it the stratospheric contribution is zero and the
     provenance says "tropospheric only". Its moments are interpolated linearly between monthly midpoints (model
     date; the file covers a run of months and a date outside them is refused) and in latitude (clamped at the end
     latitudes), then summed over model layers by height overlap (conservative; zero outside 5-40 km).
Each radiation solve: optics of A (in every source cell and month the MACv2 layer moments tau f, tau f ssa and
tau f ssa g are formed first, then interpolated linearly between monthly midpoints and bilinearly from the 1 degree grid
to each column, and distributed over model layers by height overlap of the MACv2 layers, clipped at the model surface;
Codex 14:46) plus B
(evaluated on the model layer heights) form the tropospheric optics; C forms the stratospheric optics; the
instantaneous local tropopause (WMO lapse-rate rule on model layers) assigns them: tropospheric below the tropopause
layer, stratospheric above, and in the layer holding the tropopause the ETH CMIP6 procedure, the arithmetic mean of
extinction, ssa and g of both sides, or 50 % of the extinction keeping the other quantities where one side is absent.
Optics within a component combine by optical-depth weighting (ssa) and scattering weighting (g).
=#

const CMIP_AEROSOL_BACKGROUND_SCHEMA = "readiesm_macv2_background_v1"
# The sealed MACv2 -> RRTMGP band table (band_map_macv2_rrtmgp_v1.json) and the band order it implies (Paul 14:5x).
const CMIP_AEROSOL_BAND_MAP_SHA256 = "4af3094c72101d25458a14742959dc3e953aee4ad44a823ca6535c37ae201383"
const CMIP_AEROSOL_SW_WAVELENGTHS_NM = [8000, 3470, 2790, 2320, 2050, 1780, 1460, 1270, 1000, 700, 550, 400, 300, 230]
const CMIP_AEROSOL_LW_WAVELENGTHS_NM =
    [90000, 30000, 18000, 15100, 13200, 11200, 9750, 8850, 7850, 7000, 6150, 5200, 4600, 4200, 4000, 3600]
const CMIP_AEROSOL_LABEL_TROPOSPHERIC =
    "prescribed aerosol: MACv2 pre-industrial background + MACv2-SP-1-0 (CMIP6); tropospheric only (no stratospheric file)"
# The label fixed by design v4/v5 (prescribed_aerosol_species_claude_v1/DESIGN_v5.md); never "CMIP6 prescribed aerosol".
const CMIP_AEROSOL_LABEL_FULL =
    "prescribed aerosol: MACv2 pre-industrial background + MACv2-SP-1-0 (CMIP6) + " *
    "UOEXETER-CMIP-2-2-1 (CMIP7) stratospheric optics"
const CMIP_AEROSOL_STRAT_SCHEMA = "readiesm_strat_aerosol_bands_v1"
const CMIP_AEROSOL_TROPOPAUSE_RULE =
    "WMO lapse rate on model layers (first layer above 500 hPa with lapse rate <= 2 K/km held over 2 km); ETH CMIP6 " *
    "strat-trop rule: tropospheric optics below, stratospheric above, arithmetic mean in the tropopause layer"

"""The prepared MACv2 background (host arrays in Julia order: (lon, lat, month, band) and (lon, lat, layer, month))."""
struct MACv2Background
    path::String
    sha256::String
    lon::Vector{Float64}
    lat::Vector{Float64}
    sw_wavelength_nm::Vector{Int}
    lw_wavelength_nm::Vector{Int}
    tau_fine_sw::Array{Float32, 4}
    ssa_fine_sw::Array{Float32, 4}
    g_fine_sw::Array{Float32, 4}
    tau_coarse_sw::Array{Float32, 4}
    ssa_coarse_sw::Array{Float32, 4}
    g_coarse_sw::Array{Float32, 4}
    tau_lw::Array{Float32, 4}
    ssa_lw::Array{Float32, 4}
    g_lw::Array{Float32, 4}
    frac_fine::Array{Float32, 4}
    frac_coarse::Array{Float32, 4}
    layer_bottom_m::Vector{Float64}
    layer_top_m::Vector{Float64}
end

function MACv2Background(path::AbstractString)
    isfile(path) || throw(ArgumentError("MACv2 background file not found: $path"))
    b = NCDatasets.NCDataset(path) do ds
        get(ds.attrib, "schema", "") == CMIP_AEROSOL_BACKGROUND_SCHEMA ||
            throw(ArgumentError("$path does not have schema $CMIP_AEROSOL_BACKGROUND_SCHEMA"))
        occursin(CMIP_AEROSOL_BAND_MAP_SHA256, get(ds.attrib, "band_map", "")) ||
            throw(ArgumentError("$path was not prepared with band map $CMIP_AEROSOL_BAND_MAP_SHA256"))
        r(n) = Array(ds[n].var[:, :, :, :])
        for n in ("tau_fine_sw", "ssa_fine_sw", "g_fine_sw", "tau_coarse_sw", "ssa_coarse_sw", "g_coarse_sw",
                  "tau_lw", "ssa_lw", "g_lw", "layer_fraction_fine", "layer_fraction_coarse")
            any(k -> haskey(ds[n].attrib, k), ("_FillValue", "missing_value")) &&
                throw(ArgumentError("$n in $path declares a fill value"))
        end
        (lon = Float64.(Array(ds["lon"].var[:])), lat = Float64.(Array(ds["lat"].var[:])),
         sw = Int.(Array(ds["sw_wavelength_nm"].var[:])), lw = Int.(Array(ds["lw_wavelength_nm"].var[:])),
         tfs = r("tau_fine_sw"), sfs = r("ssa_fine_sw"), gfs = r("g_fine_sw"), tcs = r("tau_coarse_sw"),
         scs = r("ssa_coarse_sw"), gcs = r("g_coarse_sw"), tl = r("tau_lw"), sl = r("ssa_lw"), gl = r("g_lw"),
         ff = r("layer_fraction_fine"), fc = r("layer_fraction_coarse"),
         zb = Float64.(Array(ds["layer_bottom_m"].var[:])), zt = Float64.(Array(ds["layer_top_m"].var[:])))
    end
    all(diff(b.lat) .> 0) && all(diff(b.lon) .> 0) || throw(ArgumentError("$path coordinates are not increasing"))
    b.sw == CMIP_AEROSOL_SW_WAVELENGTHS_NM && b.lw == CMIP_AEROSOL_LW_WAVELENGTHS_NM ||
        throw(ArgumentError("$path band wavelengths $(b.sw) / $(b.lw) differ from the sealed band order"))
    for (n, x) in (("tau", (b.tfs, b.tcs, b.tl)), ("ssa", (b.sfs, b.scs, b.sl)), ("g", (b.gfs, b.gcs, b.gl)))
        for a in x
            all(isfinite, a) || throw(ArgumentError("$path $n has nonfinite values"))
            ok = n == "tau" ? all(>=(0), a) : n == "ssa" ? all(v -> 0 <= v <= 1, a) : all(v -> 0 <= v < 1, a)
            ok || throw(ArgumentError("$path $n is outside its physical range"))
        end
    end
    for f in (b.ff, b.fc)
        all(v -> 0 <= v <= 1, f) && maximum(abs.(sum(f; dims = 3) .- 1)) <= 1.0e-4 ||
            throw(ArgumentError("$path layer fractions do not sum to 1"))
    end
    return MACv2Background(String(path), bytes2hex(open(SHA.sha256, path)), b.lon, b.lat, b.sw, b.lw, b.tfs, b.sfs,
                           b.gfs, b.tcs, b.scs, b.gcs, b.tl, b.sl, b.gl, b.ff, b.fc, b.zb, b.zt)
end

"""
Prepared stratospheric aerosol optics (component C): band moments per metre, E = extinction, S = extinction x ssa,
G = extinction x ssa x g, host arrays (band, height, lat, month) on the source latitudes, height layers and months.
"""
struct StratosphericAerosol
    path::String
    sha256::String
    lat::Vector{Float64}
    height_bottom_m::Vector{Float64}
    height_top_m::Vector{Float64}
    month_mid::Vector{DateTime}
    E_sw::Array{Float64, 4}
    S_sw::Array{Float64, 4}
    G_sw::Array{Float64, 4}
    E_lw::Array{Float64, 4}
    S_lw::Array{Float64, 4}
    G_lw::Array{Float64, 4}
end

"""A NetCDF variable as a host array with its axes in the order `names` (by dimension name, not file order)."""
function _read_by_dim_names(ds, name, names)
    v = ds[name]
    dims = Tuple(NCDatasets.dimnames(v))
    Set(dims) == Set(names) && length(dims) == length(names) ||
        throw(ArgumentError("$name has dimensions $dims, expected $names"))
    perm = Tuple(findfirst(==(n), dims) for n in names)
    return permutedims(Array(v.var[ntuple(_ -> Colon(), length(dims))...]), perm)
end

function StratosphericAerosol(path::AbstractString)
    isfile(path) || throw(ArgumentError("stratospheric aerosol file not found: $path"))
    s = NCDatasets.NCDataset(path) do ds
        get(ds.attrib, "schema", "") == CMIP_AEROSOL_STRAT_SCHEMA ||
            throw(ArgumentError("$path does not have schema $CMIP_AEROSOL_STRAT_SCHEMA"))
        tb = ds["time_bnds"]
        get(tb.attrib, "units", "") == "days since 1850-01-01" &&
            get(tb.attrib, "calendar", "") == "proleptic_gregorian" ||
            throw(ArgumentError("$path time_bnds must be days since 1850-01-01, proleptic_gregorian"))
        days = _read_by_dim_names(ds, "time_bnds", ("bnds", "time"))
        mid = [DateTime(1850, 1, 1) + Millisecond(round(Int64, (days[1, i] + days[2, i]) / 2 * 86_400_000))
               for i in axes(days, 2)]
        hb = _read_by_dim_names(ds, "height_bnds", ("bnds", "height"))
        moments = Dict{String, Array{Float64, 4}}()
        for spec in ("sw", "lw"), m in ("E", "S", "G")
            n = "$(m)_$spec"
            any(k -> haskey(ds[n].attrib, k), ("_FillValue", "missing_value")) &&
                throw(ArgumentError("$n in $path declares a fill value"))
            moments[n] = Float64.(_read_by_dim_names(ds, n, ("$(spec)_band", "height", "lat", "time")))
        end
        (; lat = Float64.(Array(ds["lat"].var[:])), hb = Float64.(hb), mid, moments)
    end
    all(diff(s.lat) .> 0) || throw(ArgumentError("$path latitudes are not increasing"))
    all(diff(s.mid) .> Day(0)) || throw(ArgumentError("$path months are not increasing"))
    bottom, top = s.hb[1, :], s.hb[2, :]
    all(top .> bottom) && all(diff(bottom) .> 0) || throw(ArgumentError("$path height bounds are not increasing"))
    for spec in ("sw", "lw")
        E, S, G = s.moments["E_$spec"], s.moments["S_$spec"], s.moments["G_$spec"]
        all(isfinite, E) && all(isfinite, S) && all(isfinite, G) ||
            throw(ArgumentError("$path $spec moments have nonfinite values"))
        # Admissible moments, as the band-optics admission (0 <= ssa <= 1, 0 <= g < 1): 0 <= S <= E (rounding slack
        # 1e-12) and 0 <= G < S where S > 0, G = 0 where S = 0 (Paul's review of ac555fa).
        all(E .>= 0) && all(S .>= 0) && all(S .<= E .* (1 + 1.0e-12)) && all(G .>= 0) &&
            all(((s, g),) -> s > 0 ? g < s : g == 0, zip(S, G)) ||
            throw(ArgumentError("$path $spec moments are outside 0 <= S <= E, 0 <= G < S (G = 0 where S = 0)"))
    end
    m = s.moments
    return StratosphericAerosol(String(path), bytes2hex(open(SHA.sha256, path)), s.lat, bottom, top, s.mid,
                                m["E_sw"], m["S_sw"], m["G_sw"], m["E_lw"], m["S_lw"], m["G_lw"])
end

"""Linear weights between the stratospheric file's monthly midpoints for `time`; a time outside them is refused."""
function _strat_time_weights(st::StratosphericAerosol, time::DateTime)
    st.month_mid[1] <= time <= st.month_mid[end] || throw(ArgumentError(
        "model time $time is outside the stratospheric aerosol months ($(st.month_mid[1]) to $(st.month_mid[end]))"))
    i0 = clamp(searchsortedlast(st.month_mid, time), 1, length(st.month_mid) - 1)
    w = Dates.value(time - st.month_mid[i0]) / Dates.value(st.month_mid[i0 + 1] - st.month_mid[i0])
    return i0, i0 + 1, w
end

"""Linear latitude weights on the stratospheric file's latitudes, clamped at the end latitudes."""
function _strat_lat_weights(glat, lat)
    y = clamp(lat, glat[1], glat[end])
    j0 = clamp(searchsortedlast(glat, y), 1, length(glat) - 1)
    return j0, j0 + 1, (y - glat[j0]) / (glat[j0 + 1] - glat[j0])
end

"""
Stratospheric band moments of one column on the model layers, `M[k, band, layer]` with k = 1, 2, 3 for E, S, G as
optical depths: the moments per metre interpolated in time and latitude, times the overlap length of each source
height layer with each model layer (`zbc`, layer bounds bottom to top).
"""
function _strat_column_moments!(M, E, S, G, st::StratosphericAerosol, lw, tw, zbc)
    j0, j1, wy = lw
    t0, t1, wt = tw
    fill!(M, 0.0)
    nlay = length(zbc) - 1
    for h in eachindex(st.height_bottom_m), r in 1:nlay
        len = min(st.height_top_m[h], zbc[r + 1]) - max(st.height_bottom_m[h], zbc[r])
        len > 0 || continue
        for b in axes(M, 2)
            for (k, A) in ((1, E), (2, S), (3, G))
                v = (1 - wt) * ((1 - wy) * A[b, h, j0, t0] + wy * A[b, h, j1, t0]) +
                    wt * ((1 - wy) * A[b, h, j0, t1] + wy * A[b, h, j1, t1])
                M[k, b, r] += v * len
            end
        end
    end
    return M
end

"""Prescribed CMIP aerosol for one model grid: sources, per-column interpolation weights and the write counter."""
mutable struct CMIPAerosol
    background::MACv2Background
    plumes::MACv2SPParameters
    fixed_year::Int
    lon::Vector{Float64}                 # model column longitude (degrees east, as SpeedyWeather)
    lat::Vector{Float64}
    hweights::Vector{NTuple{6, Float64}} # (i0, i1, j0, j1, wx, wy) on the background grid, per column
    writes::Int                          # validated writes through `_set_band_optics!`
    tropopause_layer::Vector{Int}        # at the last write, RRTMGP layer index (bottom to top) per column
    stratospheric::Union{Nothing, StratosphericAerosol}
    # Nothing without a stratospheric file, so the no-C scheme adds no array to the state (Paul, job 65731200).
    strat_lat_weights::Union{Nothing, Vector{Tuple{Int, Int, Float64}}}
end

function CMIPAerosol(background_path, plumes_path, spectral_grid::SpectralGrid; fixed_year::Integer = 0,
                     stratospheric_path::AbstractString = "")
    fixed_year >= 0 || throw(ArgumentError("aerosol_cmip_year must be zero or a positive year"))
    isempty(stratospheric_path) || fixed_year == 0 || throw(ArgumentError(
        "the stratospheric aerosol file follows the model date; it cannot be combined with aerosol_cmip_year"))
    bg = MACv2Background(background_path)
    sp = MACv2SPParameters(plumes_path)
    lond, latd = RG.get_londlatds(spectral_grid.grid)
    lon, lat = Float64.(lond), Float64.(latd)
    hw = [_bilinear_weights(bg.lon, bg.lat, lon[c], lat[c]) for c in eachindex(lon)]
    st = isempty(stratospheric_path) ? nothing : StratosphericAerosol(stratospheric_path)
    sw = isnothing(st) ? nothing : [_strat_lat_weights(st.lat, lat[c]) for c in eachindex(lat)]
    return CMIPAerosol(bg, sp, Int(fixed_year), lon, lat, hw, 0, zeros(Int, length(lon)), st, sw)
end

"""Bilinear weights on a periodic longitude axis and a clamped latitude axis."""
function _bilinear_weights(glon, glat, lon, lat)
    x = mod(lon - glon[1], 360.0) + glon[1]
    n = length(glon)
    i0 = clamp(searchsortedlast(glon, x), 1, n)
    i1 = i0 == n ? 1 : i0 + 1
    span = i0 == n ? glon[1] + 360 - glon[n] : glon[i1] - glon[i0]
    wx = x >= glon[i0] ? (x - glon[i0]) / span : (x + 360 - glon[i0]) / span
    y = clamp(lat, glat[1], glat[end])
    j0 = clamp(searchsortedlast(glat, y), 1, length(glat) - 1)
    wy = (y - glat[j0]) / (glat[j0 + 1] - glat[j0])
    return (Float64(i0), Float64(i1), Float64(j0), Float64(j0 + 1), wx, wy)
end

@inline function _bilin(a, w, m, band)
    i0, i1, j0, j1 = Int(w[1]), Int(w[2]), Int(w[3]), Int(w[4])
    wx, wy = w[5], w[6]
    return (1 - wy) * ((1 - wx) * a[i0, j0, m, band] + wx * a[i1, j0, m, band]) +
           wy * ((1 - wx) * a[i0, j1, m, band] + wx * a[i1, j1, m, band])
end

"""Monthly-midpoint weights (cyclic) for decimal-year fraction `f`: months (m0, m1) and weight of m1."""
function _month_weights(f)
    x = 12 * f - 0.5                 # month midpoints at (m - 0.5) / 12
    m0 = floor(Int, x)
    w = x - m0
    return mod(m0, 12) + 1, mod(m0 + 1, 12) + 1, w
end

"""
Instantaneous local tropopause per column (WMO lapse-rate rule on model layers, Codex 14:46): the lowest layer `r`
above 500 hPa from which the average lapse rate to EVERY higher layer within 2 km is <= 2 K/km. Coverage rule for
coarse grids: the next higher layer is always tested, even when it is more than 2 km above `r`, so a candidate is never
accepted without at least one tested interval. Returns the RRTMGP layer index (bottom to top) of the layer holding the
tropopause, or the top layer if none qualifies.
"""
function _tropopause_layers(t_lay, p_lay, z)
    nlay, ncol = size(t_lay)
    out = fill(nlay, ncol)
    for c in 1:ncol
        for r in 1:(nlay - 1)
            p_lay[r, c] < 5.0e4 || continue
            Γ(r2) = -(t_lay[r2, c] - t_lay[r, c]) / max((z[r2, c] - z[r, c]) / 1000, 1.0e-6)
            ok = Γ(r + 1) <= 2
            rr = r + 2
            while ok && rr <= nlay && z[rr, c] <= z[r, c] + 2000
                ok = Γ(rr) <= 2
                rr += 1
            end
            ok || continue
            out[c] = r
            break # Leave this column only; the next column still needs its own tropopause.
        end
    end
    return out
end

"""
MACv2 layer moments at one column for band `b`: in each of the four source cells and the two bracketing months the
layer optical depth e = tau f (f the layer fraction), its scattering e ssa and e ssa g are formed first, then weighted
bilinearly in space and linearly in time (Codex 14:46: interpolating tau and f separately and multiplying is not the
interpolated layer optical depth). Fills E, S, A (one value per MACv2 layer).
"""
function _background_layer_moments!(E, S, A, tau, ssa, g, frac, w, m0, m1, wm, b)
    fill!(E, 0.0); fill!(S, 0.0); fill!(A, 0.0)
    i0, i1, j0, j1 = Int(w[1]), Int(w[2]), Int(w[3]), Int(w[4])
    wx, wy = w[5], w[6]
    corners = ((i0, j0, (1 - wx) * (1 - wy)), (i1, j0, wx * (1 - wy)), (i0, j1, (1 - wx) * wy), (i1, j1, wx * wy))
    for (m, wt) in ((m0, 1 - wm), (m1, wm)), (i, j, wc) in corners
        ww = wt * wc
        ww == 0 && continue
        τ, ω, gg = Float64(tau[i, j, m, b]), Float64(ssa[i, j, m, b]), Float64(g[i, j, m, b])
        for l in eachindex(E)
            e = τ * frac[i, j, l, m]
            E[l] += ww * e
            S[l] += ww * e * ω
            A[l] += ww * e * ω * gg
        end
    end
    return nothing
end

"""
Optics (tau, ssa, g) from summed moments: tau = E, ssa = S / E, g = A / S; zero extinction gives (0, 1, 0). A
nonfinite or negative moment gives NaN optics, so `_set_band_optics!` refuses it instead of a comparison silently
turning it into zero (job 65679738: `NaN > 0` is false).
"""
@inline function _from_moments(E, S, A)
    (isfinite(E) && isfinite(S) && isfinite(A) && E >= 0 && S >= 0 && A >= 0) || return (NaN, NaN, NaN)
    return E > 0 ? (E, S / E, S > 0 ? A / S : 0.0) : (0.0, 1.0, 0.0)
end

"""Optics (tau, ssa, g) inside the admissible domain: finite, tau >= 0, 0 <= ssa <= 1, 0 <= g < 1."""
@inline _admissible_optics(p) = all(isfinite, p) && p[1] >= 0 && 0 <= p[2] <= 1 && 0 <= p[3] < 1

"""
Combine optics (tau, ssa, g) tuples by optical-depth and scattering weighting; any part outside the admissible domain
(nonfinite, negative tau, ssa outside [0, 1], g outside [0, 1)) gives NaN optics (Paul 15:3x).
"""
@inline function _combine_optics(parts)
    all(_admissible_optics, parts) || return (NaN, NaN, NaN)
    τ = sum(p[1] for p in parts)
    τ > 0 || return (0.0, 1.0, 0.0)
    sca = sum(p[1] * p[2] for p in parts)
    asy = sum(p[1] * p[2] * p[3] for p in parts)
    return (τ, sca / τ, sca > 0 ? asy / sca : 0.0)
end

"""
ETH CMIP6 tropopause-layer rule: arithmetic means, or 50 % of the extinction where one side is absent. Input outside
the admissible domain on either side gives NaN optics.
"""
@inline function _tropopause_mean(trop, strat)
    _admissible_optics(trop) && _admissible_optics(strat) || return (NaN, NaN, NaN)
    trop[1] > 0 && strat[1] > 0 &&
        return ((trop[1] + strat[1]) / 2, (trop[2] + strat[2]) / 2, (trop[3] + strat[3]) / 2)
    trop[1] > 0 && return (trop[1] / 2, trop[2], trop[3])
    strat[1] > 0 && return (strat[1] / 2, strat[2], strat[3])
    return (0.0, 1.0, 0.0)
end

"""
Model layer bounds in height above sea level, RRTMGP order (bottom to top), from the full-level geopotential
(SpeedyWeather order, top to bottom) and the orography: interfaces midway between full levels, the bottom at the
surface and the top extrapolated by half a layer.
"""
function _layer_heights(geopotential, orography, gravity, nlay)
    npoints, nlev = size(geopotential)
    zc = zeros(nlay, npoints); zb = zeros(nlay + 1, npoints)
    for c in 1:npoints
        for r in 1:nlev
            zc[r, c] = geopotential[c, nlev - r + 1] / gravity
        end
        for r in (nlev + 1):nlay
            zc[r, c] = zc[nlev, c] + (r - nlev) * (zc[nlev, c] - zc[nlev - 1, c])
        end
        zb[1, c] = orography[c]
        for r in 1:(nlay - 1)
            zb[r + 1, c] = (zc[r, c] + zc[r + 1, c]) / 2
        end
        zb[nlay + 1, c] = zc[nlay, c] + (zc[nlay, c] - zb[nlay, c])
    end
    return zc, zb
end

"""
Prescribed optics on host arrays `(band, layer, column)` for the current state: the tropospheric part (A + B), the
stratospheric part (C, zero without its file) and the tropopause rule.
"""
function _cmip_aerosol_optics(a::CMIPAerosol, time, zc, zb, orography, t_lay, p_lay, nbnd_sw, nbnd_lw)
    bg = a.background
    nlay, ncol = size(zc)
    length(bg.sw_wavelength_nm) == nbnd_sw && length(bg.lw_wavelength_nm) == nbnd_lw ||
        throw(ArgumentError("MACv2 background bands differ from the RRTMGP bands"))
    x = _decimal_year(time)
    frac = x - floor(x)
    year_fr = a.fixed_year > 0 ? a.fixed_year + frac : x
    m0, m1, wm = _month_weights(frac)
    tau_sw = zeros(nbnd_sw, nlay, ncol); ssa_sw = ones(nbnd_sw, nlay, ncol); g_sw = zeros(nbnd_sw, nlay, ncol)
    tau_lw = zeros(nbnd_lw, nlay, ncol); ssa_lw = ones(nbnd_lw, nlay, ncol); g_lw = zeros(nbnd_lw, nlay, ncol)
    dz = zb[2:end, :] .- zb[1:(end - 1), :]
    sp_tau, sp_ssa, sp_g = macv2sp_band_optics(a.plumes, bg.sw_wavelength_nm, orography, a.lon, a.lat, year_fr,
                                               permutedims(zc), permutedims(dz))
    nl = length(bg.layer_bottom_m)
    st = a.stratospheric
    if !isnothing(st)
        size(st.E_sw, 1) == nbnd_sw && size(st.E_lw, 1) == nbnd_lw ||
            throw(ArgumentError("stratospheric aerosol bands differ from the RRTMGP bands"))
        tw = _strat_time_weights(st, time)
    end
    Msw, Mlw = zeros(3, nbnd_sw, nlay), zeros(3, nbnd_lw, nlay)
    tp = _tropopause_layers(t_lay, p_lay, zc)
    a.tropopause_layer .= tp
    over = zeros(nl, nlay)
    Ef, Sf, Af, Ec, Sc, Ac = (zeros(nl) for _ in 1:6)
    for c in 1:ncol
        w = a.hweights[c]
        if !isnothing(st)
            zbc = view(zb, :, c)
            _strat_column_moments!(Msw, st.E_sw, st.S_sw, st.G_sw, st, a.strat_lat_weights[c], tw, zbc)
            _strat_column_moments!(Mlw, st.E_lw, st.S_lw, st.G_lw, st, a.strat_lat_weights[c], tw, zbc)
        end
        # MACv2 layer overlap with each model layer, clipped at the model surface (uniform within a MACv2 layer).
        for r in 1:nlay, l in 1:nl
            lo = max(bg.layer_bottom_m[l], zb[r, c], orography[c])
            hi = min(bg.layer_top_m[l], zb[r + 1, c])
            over[l, r] = max(hi - lo, 0.0) / (bg.layer_top_m[l] - bg.layer_bottom_m[l])
        end
        for b in 1:nbnd_sw
            _background_layer_moments!(Ef, Sf, Af, bg.tau_fine_sw, bg.ssa_fine_sw, bg.g_fine_sw, bg.frac_fine, w,
                                       m0, m1, wm, b)
            _background_layer_moments!(Ec, Sc, Ac, bg.tau_coarse_sw, bg.ssa_coarse_sw, bg.g_coarse_sw,
                                       bg.frac_coarse, w, m0, m1, wm, b)
            for r in 1:nlay
                E = S = Asum = 0.0
                for l in 1:nl
                    o = over[l, r]
                    o == 0 && continue
                    E += (Ef[l] + Ec[l]) * o; S += (Sf[l] + Sc[l]) * o; Asum += (Af[l] + Ac[l]) * o
                end
                τsp = sp_tau[b, c, r]
                E += τsp; S += τsp * sp_ssa[b, c, r]; Asum += τsp * sp_ssa[b, c, r] * sp_g[b, c, r]
                trop = _from_moments(E, S, Asum)
                strat = _from_moments(Msw[1, b, r], Msw[2, b, r], Msw[3, b, r])   # component C (zero without its file)
                o = r < tp[c] ? trop : r > tp[c] ? strat : _tropopause_mean(trop, strat)
                tau_sw[b, r, c], ssa_sw[b, r, c], g_sw[b, r, c] = o
            end
        end
        for b in 1:nbnd_lw
            # LW: the MACv2 IR total on the coarse profile (design v3).
            _background_layer_moments!(Ec, Sc, Ac, bg.tau_lw, bg.ssa_lw, bg.g_lw, bg.frac_coarse, w, m0, m1, wm, b)
            for r in 1:nlay
                E = sum(Ec[l] * over[l, r] for l in 1:nl)
                S = sum(Sc[l] * over[l, r] for l in 1:nl)
                Asum = sum(Ac[l] * over[l, r] for l in 1:nl)
                trop = _from_moments(E, S, Asum)
                strat = _from_moments(Mlw[1, b, r], Mlw[2, b, r], Mlw[3, b, r])
                o = r < tp[c] ? trop : r > tp[c] ? strat : _tropopause_mean(trop, strat)
                tau_lw[b, r, c], ssa_lw[b, r, c], g_lw[b, r, c] = o
            end
        end
    end
    return (; tau_sw, ssa_sw, g_sw, tau_lw, ssa_lw, g_lw)
end

"""
The only PA2 write path for band optics (Paul, 14:xx): checks the host optics (shape, finite, tau >= 0, 0 <= ssa <= 1,
0 <= g < 1) and the zero-MERRA invariant on the solver's own architecture, then copies them into the solver and counts
the write. Throws before anything is copied if a check fails.
"""
function _set_band_optics!(a::CMIPAerosol, solver, optics)
    bo = solver.as.aerosol_state.band_optics
    isnothing(bo) && throw(ArgumentError("the RRTMGP solver has no band-optics arrays"))
    keys(optics) == keys(bo) || throw(ArgumentError("band optics fields $(keys(optics)) differ from $(keys(bo))"))
    # Convert to the destination precision first and check what will be stored (Paul 14:5x: a finite Float64 can
    # become nonfinite in Float32).
    optics = (; (name => eltype(getproperty(bo, name)).(getproperty(optics, name)) for name in keys(bo))...)
    for name in keys(bo)
        x = getproperty(optics, name)
        size(x) == size(getproperty(bo, name)) ||
            throw(ArgumentError("band optics $name has size $(size(x)), expected $(size(getproperty(bo, name)))"))
        all(isfinite, x) || throw(ArgumentError("band optics $name has a nonfinite value"))
        s = String(name)
        ok = startswith(s, "tau") ? all(>=(0), x) : startswith(s, "ssa") ? all(v -> 0 <= v <= 1, x) :
             all(v -> 0 <= v < 1, x)
        ok || throw(ArgumentError("band optics $name is outside its physical range"))
    end
    all(iszero, solver.as.aerosol_state.aero_mass) || throw(ArgumentError(
        "MERRA aerosol mass must be zero when prescribed band optics own the aerosol (no double counting)"))
    for name in keys(bo)
        dst = getproperty(bo, name)
        copyto!(dst, getproperty(optics, name))
    end
    a.writes += 1
    return nothing
end

"""Compute and write the prescribed optics before a flux solve. A no-op without the CMIP aerosol scheme."""
_apply_cmip_aerosol!(::Nothing, solver, vars, model; upper_cap = nothing) = nothing

function _apply_cmip_aerosol!(a::CMIPAerosol, solver, vars, model; upper_cap = nothing)
    nlay = solver.grid_params.nlay
    geopotential = Array(vars.grid.geopotential.data)
    orography = Float64.(vec(Array(model.orography.orography.data)))
    zc, zb = _upper_cap_layer_heights(
        upper_cap, Float64.(geopotential), orography, model.planet.gravity,
        solver, model.geometry.nlayers,
    )
    t_lay = Float64.(Array(RRTMGP.layer_temperature(solver)))
    p_lay = Float64.(Array(RRTMGP.layer_pressure(solver)))
    size(t_lay, 1) == nlay || (t_lay = vcat(t_lay, repeat(t_lay[end:end, :], nlay - size(t_lay, 1))))
    size(p_lay, 1) == nlay || (p_lay = vcat(p_lay, repeat(p_lay[end:end, :], nlay - size(p_lay, 1))))
    optics = _cmip_aerosol_optics(a, vars.prognostic.clock.time, zc, zb, orography, t_lay, p_lay,
                                  solver.lookups.nbnd_sw, solver.lookups.nbnd_lw)
    if !isnothing(upper_cap)
        # First cap contract: no aerosol optics above the 5 hPa split. Preserve
        # the hydrostatic geometry for the dynamic top layer, but do not let
        # extrapolated stratospheric optics populate the two provisional caps.
        first_cap = model.geometry.nlayers + 1
        for name in keys(optics)
            values = getproperty(optics, name)
            fill!(@view(values[:, first_cap:end, :]), startswith(String(name), "ssa") ? 1.0 : 0.0)
        end
    end
    _set_band_optics!(a, solver, optics)
    return nothing
end

"""The scheme for `config`, or `nothing` unless forcing.aerosol_scheme = :cmip_prescribed."""
_cmip_aerosol(config, spectral_grid) = config.forcing.aerosol_scheme == :cmip_prescribed ?
    CMIPAerosol(config.forcing.aerosol_cmip_background_path, config.forcing.aerosol_cmip_plumes_path, spectral_grid;
                fixed_year = config.forcing.aerosol_cmip_year,
                stratospheric_path = config.forcing.aerosol_cmip_stratospheric_path) : nothing

"""Restart identity of the scheme: source hashes and year mode."""
_cmip_aerosol_identity(::Nothing) = nothing
_cmip_aerosol_identity(a::CMIPAerosol) =
    (background_sha256 = a.background.sha256, plumes_sha256 = a.plumes.sha256, fixed_year = a.fixed_year,
     stratospheric_sha256 = isnothing(a.stratospheric) ? "" : a.stratospheric.sha256)

"""Reject a checkpoint whose CMIP aerosol identity differs from the reconstructed radiation."""
function _validate_cmip_aerosol_restart(radiation, state)
    saved = hasproperty(state, :cmip_aerosol) ? state.cmip_aerosol : nothing
    current = _cmip_aerosol_identity(radiation.cmip_aerosol)
    saved == current || throw(ArgumentError(
        "checkpoint and reconstructed radiation use different CMIP aerosol (checkpoint $(something(saved, "none")), " *
        "configuration $(something(current, "none")))"))
    return nothing
end

"""Diagnostics provenance of the CMIP aerosol scheme, written only when it is requested."""
function _cmip_aerosol_provenance(radiation, config, rrtmgp_active)
    config.forcing.aerosol_scheme == :cmip_prescribed || return (;)
    a = rrtmgp_active ? radiation.cmip_aerosol : nothing
    isnothing(a) && return (cmip_aerosol = "requested_but_not_active",)
    st = a.stratospheric
    return (cmip_aerosol = isnothing(st) ? CMIP_AEROSOL_LABEL_TROPOSPHERIC : CMIP_AEROSOL_LABEL_FULL,
            cmip_aerosol_stratospheric_path = isnothing(st) ? "" : st.path,
            cmip_aerosol_stratospheric_sha256 = isnothing(st) ? "" : st.sha256,
            cmip_aerosol_background_path = a.background.path, cmip_aerosol_background_sha256 = a.background.sha256,
            cmip_aerosol_plumes_path = a.plumes.path, cmip_aerosol_plumes_sha256 = a.plumes.sha256,
            cmip_aerosol_year_mode = a.fixed_year > 0 ? "fixed_year_$(a.fixed_year)" : "model_date",
            cmip_aerosol_band_rule = "band-representative wavelengths (MACv2 convention), band map 4af3094c",
            cmip_aerosol_tropopause_rule = CMIP_AEROSOL_TROPOPAUSE_RULE,
            cmip_aerosol_merra_mass = "zero, checked at every write",
            cmip_aerosol_validated_writes = a.writes)
end
