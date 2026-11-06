#=
MACv2-SP anthropogenic aerosol optical properties (PA2 component B), a Julia port of the CMIP6 reference routine
`sp_aop_profile` in mo_simple_plumes_v1.f90 (input4MIPs MPI-M-MACv2-SP-1-0, v20170201, sha256 1fc26c89; Stevens et
al. 2017, GMD 10, 433-452, doi:10.5194/gmd-10-433-2017). The parameters come from MACv2.0-SP_v1.nc (sha256 b99f1372).

The port follows the reference line by line: nine plumes of two rotated anisotropic Gaussian features each; feature
amplitudes = year weight x weekly annual cycle x feature weight x AOD at source; a beta-function vertical profile in
eta = z / 15 km, set to zero below the orography and normalised over the column; Angstrom wavelength scaling of the
AOD; the wavelength dependence of SSA and asymmetry (lfactor = min(1, 700/lambda)); optical-depth-weighted
superposition of the plumes; and the Twomey factor dN/N from the anthropogenic and fine-mode background column AODs.
The Twomey factor is computed but not used by PA2 (a later opt-in). Following Stevens et al. (2017, sec. 2.3), the
effect at wavelengths above 3000 nm is neglected: `macv2sp_band_optics` gives zero optical depth there.
=#

const MACV2SP_NPLUMES = 9
const MACV2SP_NFEATURES = 2
const MACV2SP_NTIMES = 52
const MACV2SP_FIRST_YEAR = 1850
const MACV2SP_NYEARS = 251

"""MACv2-SP plume parameters read from MACv2.0-SP_v1.nc, as `sp_setup` reads them."""
struct MACv2SPParameters
    path::String
    sha256::String
    plume_lat::Vector{Float64}
    plume_lon::Vector{Float64}
    beta_a::Vector{Float64}
    beta_b::Vector{Float64}
    aod_spmx::Vector{Float64}
    aod_fmbg::Vector{Float64}
    asy550::Vector{Float64}
    ssa550::Vector{Float64}
    angstrom::Vector{Float64}
    sig_lon_E::Matrix{Float64}     # (feature, plume)
    sig_lon_W::Matrix{Float64}
    sig_lat_E::Matrix{Float64}
    sig_lat_W::Matrix{Float64}
    theta::Matrix{Float64}
    ftr_weight::Matrix{Float64}
    year_weight::Matrix{Float64}   # (year, plume)
    ann_cycle::Array{Float64, 3}   # (feature, week, plume)
    last_year::Int                 # last year with year weights (MACv2-SP-1-0 is the CMIP6 historical file: 2016)
end

function MACv2SPParameters(path::AbstractString)
    isfile(path) || throw(ArgumentError("MACv2-SP parameter file not found: $path"))
    p = NCDatasets.NCDataset(path) do ds
        dims = (ds.dim["plume_number"], ds.dim["plume_feature"], ds.dim["year_fr"], ds.dim["years"])
        dims == (MACV2SP_NPLUMES, MACV2SP_NFEATURES, MACV2SP_NTIMES, MACV2SP_NYEARS) ||
            throw(ArgumentError("MACv2-SP file $path is improperly dimensioned: $dims"))
        Int(ds["years"].var[1]) == MACV2SP_FIRST_YEAR ||
            throw(ArgumentError("MACv2-SP years must start in $MACV2SP_FIRST_YEAR"))
        # `.var` reads raw values; every parameter declares a _FillValue (9.96921e36), so none may equal it.
        get1(n) = Float64.(Array(ds[n].var[:]))
        get2(n) = Float64.(Array(ds[n].var[:, :]))
        out = (plume_lat = get1("plume_lat"), plume_lon = get1("plume_lon"), beta_a = get1("beta_a"),
               beta_b = get1("beta_b"), aod_spmx = get1("aod_spmx"), aod_fmbg = get1("aod_fmbg"),
               asy550 = get1("asy550"), ssa550 = get1("ssa550"), angstrom = get1("angstrom"),
               sig_lon_E = get2("sig_lon_E"), sig_lon_W = get2("sig_lon_W"), sig_lat_E = get2("sig_lat_E"),
               sig_lat_W = get2("sig_lat_W"), theta = get2("theta"), ftr_weight = get2("ftr_weight"),
               year_weight = get2("year_weight"), ann_cycle = Float64.(Array(ds["ann_cycle"].var[:, :, :])))
        for (name, v) in pairs(out)
            name == :year_weight && continue
            all(isfinite, v) && all(x -> abs(x) < 1.0e30, v) ||
                throw(ArgumentError("MACv2-SP $name in $path has a fill or nonfinite value"))
        end
        out
    end
    # The CMIP6 historical file fills year_weight after its last historical year (2017-2100 are the scenario files'
    # job). Fill is accepted only as one trailing block of years common to all plumes; those years are refused at use.
    yw = p.year_weight
    valid = [all(v -> isfinite(v) && abs(v) < 1.0e30, view(yw, iy, :)) for iy in axes(yw, 1)]
    last_valid = findlast(valid)
    !isnothing(last_valid) && all(valid[1:last_valid]) && !any(valid[(last_valid + 1):end]) ||
        throw(ArgumentError("MACv2-SP year_weight in $path has fill values that are not one trailing block of years"))
    any(iy -> any(v -> !isfinite(v) || abs(v) >= 1.0e30, view(yw, iy, :)), 1:last_valid) &&
        throw(ArgumentError("MACv2-SP year_weight in $path has fill values inside its valid years"))
    # NCDatasets reverses the NetCDF dimension order, so (nplumes, nfeatures) on disk is (nfeatures, nplumes) here,
    # as the Fortran arrays; year_weight is (years, plumes) and ann_cycle (features, weeks, plumes), as in Fortran.
    size(p.sig_lon_E) == (MACV2SP_NFEATURES, MACV2SP_NPLUMES) &&
        size(p.year_weight) == (MACV2SP_NYEARS, MACV2SP_NPLUMES) &&
        size(p.ann_cycle) == (MACV2SP_NFEATURES, MACV2SP_NTIMES, MACV2SP_NPLUMES) ||
        throw(ArgumentError("MACv2-SP arrays in $path have unexpected shapes"))
    return MACv2SPParameters(String(path), bytes2hex(open(SHA.sha256, path)), p.plume_lat, p.plume_lon, p.beta_a,
                             p.beta_b, p.aod_spmx, p.aod_fmbg, p.asy550, p.ssa550, p.angstrom, p.sig_lon_E,
                             p.sig_lon_W, p.sig_lat_E, p.sig_lat_W, p.theta, p.ftr_weight, p.year_weight, p.ann_cycle,
                             MACV2SP_FIRST_YEAR + last_valid - 1)
end

"""Feature time weights for fractional year `year_fr`, as `set_time_weight`: (anthropogenic, background)."""
function _macv2sp_time_weights(p::MACv2SPParameters, year_fr::Real)
    iyear = floor(Int, year_fr) - (MACV2SP_FIRST_YEAR - 1)
    iweek = floor(Int, (year_fr - floor(year_fr)) * MACV2SP_NTIMES) + 1
    1 <= iweek <= MACV2SP_NTIMES && 1 <= iyear <= MACV2SP_NYEARS ||
        throw(ArgumentError("MACv2-SP time $year_fr is outside 1850-2100"))
    floor(Int, year_fr) <= p.last_year ||
        throw(ArgumentError("MACv2-SP time $year_fr is after the last year with weights in $(p.path) ($(p.last_year))"))
    tw = [p.year_weight[iyear, ip] * p.ann_cycle[f, iweek, ip] for f in 1:MACV2SP_NFEATURES, ip in 1:MACV2SP_NPLUMES]
    tw_bg = [p.ann_cycle[f, iweek, ip] for f in 1:MACV2SP_NFEATURES, ip in 1:MACV2SP_NPLUMES]
    return tw, tw_bg
end

"""
    macv2sp_aop_profile(p, lambda, oro, lon, lat, year_fr, z, dz) -> (aod, ssa, asy, dNovrN)

The reference `sp_aop_profile` for columns `c` with orography `oro[c]` (m), longitude and latitude (degrees), heights
`z[c, k]` (m above sea level) and thicknesses `dz[c, k]` (m) at wavelength `lambda` (nm). Returns the layer optical
depth, single-scattering albedo and asymmetry, each `(ncol, nlev)`, and the Twomey factor per column.
"""
function macv2sp_aop_profile(p::MACv2SPParameters, lambda::Real, oro, lon, lat, year_fr::Real, z, dz)
    ncol, nlev = size(z)
    tw, tw_bg = _macv2sp_time_weights(p, year_fr)
    aod = zeros(ncol, nlev); ssa = zeros(ncol, nlev); asy = zeros(ncol, nlev)
    z_beta = [z[c, k] >= oro[c] ? 1.0 : 0.0 for c in 1:ncol, k in 1:nlev]
    eta = [max(0.0, min(1.0, z[c, k] / 15000.0)) for c in 1:ncol, k in 1:nlev]
    caod_sp = zeros(ncol); caod_bg = fill(0.02, ncol)
    prof = zeros(ncol, nlev)
    cw_an = zeros(ncol); cw_bg = zeros(ncol); ssa_c = zeros(ncol); asy_c = zeros(ncol)
    for ip in 1:MACV2SP_NPLUMES
        beta_sum = zeros(ncol)
        for k in 1:nlev, c in 1:ncol
            prof[c, k] = (eta[c, k]^(p.beta_a[ip] - 1) * (1 - eta[c, k])^(p.beta_b[ip] - 1)) * dz[c, k]
            beta_sum[c] += prof[c, k]
        end
        for k in 1:nlev, c in 1:ncol
            prof[c, k] = (prof[c, k] / beta_sum[c]) * z_beta[c, k]
        end
        for c in 1:ncol
            delta_lat = lat[c] - p.plume_lat[ip]
            delta_lon = lon[c] - p.plume_lon[ip]
            delta_lon_t = ip == 1 ? 260.0 : 180.0
            delta_lon = abs(delta_lon) > delta_lon_t ? delta_lon - copysign(360.0, delta_lon) : delta_lon
            east = delta_lon > 0
            a1 = 0.5 / (east ? p.sig_lon_E[1, ip] : p.sig_lon_W[1, ip])^2
            b1 = 0.5 / (east ? p.sig_lat_E[1, ip] : p.sig_lat_W[1, ip])^2
            a2 = 0.5 / (east ? p.sig_lon_E[2, ip] : p.sig_lon_W[2, ip])^2
            b2 = 0.5 / (east ? p.sig_lat_E[2, ip] : p.sig_lat_W[2, ip])^2
            lon1 = cos(p.theta[1, ip]) * delta_lon + sin(p.theta[1, ip]) * delta_lat
            lat1 = -sin(p.theta[1, ip]) * delta_lon + cos(p.theta[1, ip]) * delta_lat
            lon2 = cos(p.theta[2, ip]) * delta_lon + sin(p.theta[2, ip]) * delta_lat
            lat2 = -sin(p.theta[2, ip]) * delta_lon + cos(p.theta[2, ip]) * delta_lat
            e1 = exp(-(a1 * lon1^2 + b1 * lat1^2))
            e2 = exp(-(a2 * lon2^2 + b2 * lat2^2))
            f1 = tw[1, ip] * p.ftr_weight[1, ip] * e1
            f2 = tw[2, ip] * p.ftr_weight[2, ip] * e2
            f3 = tw_bg[1, ip] * p.ftr_weight[1, ip] * e1
            f4 = tw_bg[2, ip] * p.ftr_weight[2, ip] * e2
            cw_an[c] = f1 * p.aod_spmx[ip] + f2 * p.aod_spmx[ip]
            cw_bg[c] = f3 * p.aod_fmbg[ip] + f4 * p.aod_fmbg[ip]
            lf = min(1.0, 700.0 / lambda)
            ssa_c[c] = (p.ssa550[ip] * lf^4) / ((p.ssa550[ip] * lf^4) + ((1 - p.ssa550[ip]) * lf))
            asy_c[c] = p.asy550[ip] * sqrt(lf)
        end
        lfactor = exp(-p.angstrom[ip] * log(lambda / 550.0))
        for k in 1:nlev, c in 1:ncol
            aod_550 = prof[c, k] * cw_an[c]
            aod_lmd = aod_550 * lfactor
            caod_sp[c] += aod_550
            caod_bg[c] += prof[c, k] * cw_bg[c]
            asy[c, k] += aod_lmd * ssa_c[c] * asy_c[c]
            ssa[c, k] += aod_lmd * ssa_c[c]
            aod[c, k] += aod_lmd
        end
    end
    tiny = floatmin(Float32)   # the reference's TINY(1.) for its default (single) REAL
    for k in 1:nlev, c in 1:ncol
        asy[c, k] = ssa[c, k] > tiny ? asy[c, k] / ssa[c, k] : 0.0
        ssa[c, k] = aod[c, k] > tiny ? ssa[c, k] / aod[c, k] : 1.0
    end
    dNovrN = [log(1000.0 * (caod_sp[c] + caod_bg[c]) + 1.0) / log(1000.0 * caod_bg[c] + 1.0) for c in 1:ncol]
    return aod, ssa, asy, dNovrN
end

"""
    macv2sp_band_optics(p, wavelengths_nm, oro, lon, lat, year_fr, z, dz) -> (tau, ssa, g)

MACv2-SP optics at the representative wavelength of each band, `(nband, ncol, nlev)`. Bands whose wavelength exceeds
3000 nm get zero optical depth (Stevens et al. 2017, sec. 2.3), with ssa 1 and g 0.
"""
function macv2sp_band_optics(p::MACv2SPParameters, wavelengths_nm, oro, lon, lat, year_fr, z, dz)
    ncol, nlev = size(z)
    nb = length(wavelengths_nm)
    tau = zeros(nb, ncol, nlev); ssa = ones(nb, ncol, nlev); g = zeros(nb, ncol, nlev)
    for (b, wl) in enumerate(wavelengths_nm)
        wl > 3000 && continue
        a, s, y, _ = macv2sp_aop_profile(p, wl, oro, lon, lat, year_fr, z, dz)
        tau[b, :, :] .= a; ssa[b, :, :] .= s; g[b, :, :] .= y
    end
    return tau, ssa, g
end
