module OceanOverturning

# Basin overturning streamfunction and northward heat transport by the staircase method (Claude, 2026-10-29,
# DESIGN.md phase 2b). Pure host functions. Index conventions are Oceananigans': (i, j, k), k = 1 at the bottom,
# x periodic, u[i, j, k] on the face between cells i-1 and i, v[i, j, k] on the face between rows j-1 and j.

export Basins, build_basins, Line, build_lines, line_transports, streamfunction

const BASIN_NAMES = ("atlantic_arctic_ocean", "indian_pacific_ocean", "global_ocean")

"""
Basin masks on the surface tracer grid. Atlantic-Arctic and Indo-Pacific come from a flood fill over wet surface cells
north of `south_cut` (CMIP: 34 S) and south of `north_limit`, started at one seed cell each. Neighbours are the four
side neighbours, periodic in x, with no fold connectivity (cells at or north of `north_limit` are outside).
"""
struct Basins
    masks::Vector{BitMatrix}       # atlantic_arctic, indian_pacific, global
end

function flood(wet, allowed, seed)
    nx, ny = size(wet)
    mask = falses(nx, ny)
    allowed[seed...] && wet[seed...] || error("basin seed $(seed) is not an allowed wet cell")
    stack = [seed]
    mask[seed...] = true
    while !isempty(stack)
        i, j = pop!(stack)
        for (di, dj) in ((1, 0), (-1, 0), (0, 1), (0, -1))
            ii, jj = mod1(i + di, nx), j + dj
            1 <= jj <= ny || continue
            (wet[ii, jj] && allowed[ii, jj] && !mask[ii, jj]) || continue
            mask[ii, jj] = true
            push!(stack, (ii, jj))
        end
    end
    return mask
end

function build_basins(wet_surface::AbstractMatrix{Bool}, latitude; atlantic_seed, pacific_seed,
                      south_cut=-34.0, north_limit=65.0)
    size(latitude) == size(wet_surface) || throw(DimensionMismatch("basin latitude"))
    allowed = (latitude .>= south_cut) .& (latitude .< north_limit)
    atlantic = flood(wet_surface, allowed, atlantic_seed)
    pacific = flood(wet_surface, allowed, pacific_seed)
    any(atlantic .& pacific) && error("Atlantic and Indo-Pacific basins connect north of $(south_cut)")
    return Basins([atlantic, pacific, BitMatrix(wet_surface)])
end

"""One staircase line: the signed faces between basin cells south of `latitude` and basin cells north of it."""
struct Line
    latitude::Float64
    basin::Int
    xfaces::Vector{NTuple{4,Int}}  # (i_face, j, i_a, i_b): the face u[i_face, j] between cells i_a and i_b
    xsign::Vector{Int}
    yfaces::Vector{NTuple{2,Int}}  # (i, j_face), sign +1 if flux j -> j+1 is northward
    ysign::Vector{Int}
end

"""
Build the lines for every target latitude and basin. A cell is south when its centre latitude is below the target.
A face counts when both its cells are in the basin and they lie on opposite sides. Lines whose target latitude is
not below every wet cell of the top row are dropped, so no face crosses the tripolar fold.
"""
function build_lines(basins::Basins, latitude, targets)
    nx, ny = size(latitude)
    wet_top = basins.masks[3][:, ny]
    fold_limit = any(wet_top) ? minimum(latitude[wet_top, ny]) : Inf
    lines = Line[]
    for (b, mask) in enumerate(basins.masks), phi in targets
        phi < fold_limit || continue
        south = (latitude .< phi) .& mask
        north = (latitude .>= phi) .& mask
        xf, xs, yf, ys = NTuple{4,Int}[], Int[], NTuple{2,Int}[], Int[]
        for j in 1:ny, i in 1:nx
            ip = mod1(i + 1, nx)                       # face between cell i and i+1 is u[ip, j]
            if south[i, j] && north[ip, j]
                push!(xf, (ip, j, i, ip)); push!(xs, +1)
            elseif north[i, j] && south[ip, j]
                push!(xf, (ip, j, ip, i)); push!(xs, -1)
            end
            if j < ny                                  # face between rows j and j+1 is v[i, j+1]
                if south[i, j] && north[i, j + 1]
                    push!(yf, (i, j + 1)); push!(ys, +1)
                elseif north[i, j] && south[i, j + 1]
                    push!(yf, (i, j + 1)); push!(ys, -1)
                end
            end
        end
        push!(lines, Line(phi, b, xf, xs, yf, ys))
    end
    return lines, fold_limit
end

"""
Per-level northward volume transport (m3 s-1) and northward heat transport (W, with rho_cp) through each line.
`u`, `v` are face velocities, `dz` the current cell thickness at centres, `theta` the cell temperature, `wet` the 3D
wet mask, `dy_x` the length of each x face (at u[i, j]) and `dx_y` of each y face (at v[i, j]). A face's thickness is
the mean of its two cells' thickness. Its temperature is the mean of its two cells (centred, not the WENO flux).
"""
function line_transports(lines, u, v, dz, theta, wet, dy_x, dx_y, rho_cp)
    nz = size(dz, 3)
    volume = zeros(length(lines), nz)
    heat = zeros(length(lines))
    for (n, line) in enumerate(lines)
        for ((iface, j, ia, ib), s) in zip(line.xfaces, line.xsign)
            for k in 1:nz
                (wet[ia, j, k] && wet[ib, j, k]) || continue
                q = s * Float64(u[iface, j, k]) * Float64(dy_x[iface, j]) * 0.5 * (Float64(dz[ia, j, k]) + Float64(dz[ib, j, k]))
                volume[n, k] += q
                heat[n] += rho_cp * q * 0.5 * (Float64(theta[ia, j, k]) + Float64(theta[ib, j, k]))
            end
        end
        for ((i, jface), s) in zip(line.yfaces, line.ysign)
            ja, jb = jface - 1, jface
            for k in 1:nz
                (wet[i, ja, k] && wet[i, jb, k]) || continue
                q = s * Float64(v[i, jface, k]) * Float64(dx_y[i, jface]) * 0.5 * (Float64(dz[i, ja, k]) + Float64(dz[i, jb, k]))
                volume[n, k] += q
                heat[n] += rho_cp * q * 0.5 * (Float64(theta[i, ja, k]) + Float64(theta[i, jb, k]))
            end
        end
    end
    return volume, heat
end

"""Overturning mass streamfunction (kg s-1) at the top face of each level: rho_0 times the northward transport
summed from the surface down to that face, so psi = 0 at the surface (CMIP sign)."""
function streamfunction(volume, rho0)
    nlines, nz = size(volume)
    psi = zeros(nlines, nz + 1)                    # interface index 1 = bottom, nz + 1 = surface
    for n in 1:nlines
        acc = 0.0
        for k in nz:-1:1
            psi[n, k + 1] = rho0 * acc
            acc += volume[n, k]
        end
        psi[n, 1] = rho0 * acc
    end
    return psi
end

end
