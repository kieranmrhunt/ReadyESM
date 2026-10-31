# Offline east/north ocean velocity from the writer's annual native face velocities (Claude, 2026-10-29).
# Offline analysis code. It builds the ocean's horizontal tripolar grid on the CPU, puts the annual mean native u and v
# into Fields with the grid's own boundary conditions, and applies the O-015 observer (commit 019aae5) level by level:
# face-to-centre average, then rotation to east and north. Rotation is linear and the grid does not move horizontally,
# so rotating the annual mean equals the annual mean of the rotated velocity.
#
# A self-test runs first. Two non-divergent solid-body flows are built from a streamfunction at cell corners, so the
# native face velocities are exact on the C grid. The observer must return the analytic east/north velocity. The second
# flow crosses the north fold, so a wrong fold sign shows as large errors in the top rows.
#
# The self-test is judged on the active ocean cells of MASK_FILE (any writer file on this grid), without the southern
# boundary row, where the default wall condition sets the normal face velocity to zero.
#
# usage: julia --project=CANDIDATE rotate_native_velocity.jl NX NY NZ DEPTH_M OUT_DIR MASK_FILE ANNUAL_FILE...
#        (OUT_DIR must not exist; with no ANNUAL_FILE only the self-test runs)
using NCDatasets, SHA, TOML
import Oceananigans
using KernelAbstractions: @kernel, @index
const O = Oceananigans

const OBSERVER = joinpath(@__DIR__, "observers", "geographic_velocity_diagnostics.jl")
const OBSERVER_SHA256 = "070be90b4b031d4ad92d7cf66e0a9f65ed57254b58d28b640c91a1893923bca9"
bytes2hex(sha256(read(OBSERVER))) == OBSERVER_SHA256 || error("O-015 observer file changed")
include(OBSERVER)

filesha(path) = bytes2hex(open(sha256, path))

function build_grid(nx, ny, nz, depth)
    z = O.ExponentialDiscretization(nz, -Float32(depth), zero(Float32); mutable=true)
    return O.TripolarGrid(O.CPU(), Float32; size=(nx, ny, nz), z, halo=(5, 5, 4))
end

node(f, i, j, grid, lx, ly) = Float64(f(i, j, 1, grid, lx, ly, O.Center()))

# The model's velocities get their boundary conditions from regularize_field_boundary_conditions with the names :u and
# :v, which gives the north fold sign -1. A plain XFaceField(grid) gets +1, which is wrong for a vector component.
function velocity_fields(grid)
    regular(name) = O.BoundaryConditions.regularize_field_boundary_conditions(O.FieldBoundaryConditions(), grid, name)
    u = O.XFaceField(grid; boundary_conditions=regular(:u))
    v = O.YFaceField(grid; boundary_conditions=regular(:v))
    u.boundary_conditions.north.condition == -1 && v.boundary_conditions.north.condition == -1 ||
        error("velocity fold sign is not -1")
    return u, v
end

function rotate(grid, uarr, varr)
    u, v = velocity_fields(grid)
    size(O.interior(u)) == size(uarr) && size(O.interior(v)) == size(varr) ||
        throw(DimensionMismatch("native u $(size(uarr)), v $(size(varr)) against Fields $(size(O.interior(u)))"))
    O.interior(u) .= uarr
    O.interior(v) .= varr
    nx, ny, nz = size(grid)
    east, north = zeros(nx, ny, nz), zeros(nx, ny, nz)
    for k in 1:nz
        r = _geographic_surface_velocity(u, v; k)
        east[:, :, k] .= r.eastward
        north[:, :, k] .= r.northward
    end
    return east, north
end

# Streamfunction flow: u = -dpsi/dy on x faces, v = dpsi/dx on y faces, psi at (Face, Face) corners.
function streamfunction_flow(grid, psi_of)
    nx, ny, nz = size(grid)
    F, C = O.Face(), O.Center()
    psi = [psi_of(node(O.λnode, i, j, grid, F, F), node(O.φnode, i, j, grid, F, F)) for i in 1:(nx + 1), j in 1:(ny + 1)]
    uarr = zeros(Float32, nx, ny, nz)
    varr = zeros(Float32, nx, ny, nz)
    for k in 1:nz, j in 1:ny, i in 1:nx
        uarr[i, j, k] = -(psi[i, j + 1] - psi[i, j]) / O.Operators.Δyᶠᶜᶜ(i, j, k, grid)
        varr[i, j, k] = (psi[i + 1, j] - psi[i, j]) / O.Operators.Δxᶜᶠᶜ(i, j, k, grid)
    end
    return uarr, varr
end

function self_test(grid, active)
    nx, ny, _ = size(grid)
    C = O.Center()
    R, U = Float64(grid.radius), 1.0
    lon = [node(O.λnode, i, j, grid, C, C) for i in 1:nx, j in 1:ny]
    lat = [node(O.φnode, i, j, grid, C, C) for i in 1:nx, j in 1:ny]
    # Cells near the two northern grid singularities (on land in the model) have strongly distorted cells.
    mapping = grid.conformal_mapping
    poles = [(Float64(mapping.first_pole_longitude) + d, Float64(mapping.north_poles_latitude)) for d in (0.0, 180.0)]
    println("NORTH_POLES ", poles)
    distance(l1, p1, l2, p2) = acosd(clamp(sind(p1) * sind(p2) + cosd(p1) * cosd(p2) * cosd(l1 - l2), -1, 1))
    pole_distance = [minimum(distance(lon[i, j], lat[i, j], pl, pp) for (pl, pp) in poles) for i in 1:nx, j in 1:ny]
    far = pole_distance .> 5.0
    judged = active .& [j > 1 for i in 1:nx, j in 1:ny]
    worst(err, mask) = (I = argmax(ifelse.(mask .& isfinite.(err), err, -Inf));
        (i=I[1], j=I[2], lon=lon[I], lat=lat[I], pole_distance=pole_distance[I], error=err[I]))
    flows = (
        ("zonal solid body", (l, p) -> -U * R * sind(p), (l, p) -> (U * cosd(p), 0.0)),
        ("cross-polar solid body", (l, p) -> -U * R * cosd(p) * cosd(l), (l, p) -> (-U * sind(p) * cosd(l), U * sind(l))),
    )
    results = Dict{String,Any}()
    pass = true
    for (name, psi_of, exact) in flows
        uarr, varr = streamfunction_flow(grid, psi_of)
        east, north = rotate(grid, uarr, varr)
        err = [hypot(east[i, j, 1] - exact(lon[i, j], lat[i, j])[1], north[i, j, 1] - exact(lon[i, j], lat[i, j])[2])
               for i in 1:nx, j in 1:ny]
        far_err = maximum(err[far .& isfinite.(err)])
        ocean_err = maximum(err[judged])
        top_err = maximum(err[:, (ny - 2):ny][judged[:, (ny - 2):ny]])
        south_err = maximum(err[:, 1])
        println("SELF_TEST_WORST ", name, " ocean ", worst(err, judged), " any ", worst(err, trues(nx, ny)))
        # Oceananigans' rotation_angle takes the grid direction from corner latitude differences, which lose
        # accuracy close to the geographic pole (measured 29 Oct: 0.0006, 0.0023, 0.014, 0.054 and 0.00001 for a
        # speed of 1 in these bands). The limits sit above those values and catch gross faults such as a wrong sign.
        bands = [(name=b, error=maximum(err[judged .& (lo .< lat .<= hi)]; init=0.0), limit=limit)
                 for (b, lo, hi, limit) in (("south_of_80N", -90, 80, 0.002), ("80N_85N", 80, 85, 0.005),
                                            ("85N_88N", 85, 88, 0.025), ("88N_89N", 88, 89, 0.1),
                                            ("north_of_89N", 89, 90, 0.1))]
        println("SELF_TEST_BANDS ", name, " ", bands)
        south_of_cap = maximum(err[judged .& (lat .<= 88)])
        cap = maximum(err[judged .& (lat .> 88)]; init=0.0)
        levels_equal = all(isequal(east[:, :, k], east[:, :, 1]) && isequal(north[:, :, k], north[:, :, 1])
                           for k in axes(east, 3))
        nonfinite = .!isfinite.(err)
        ok = all(b.error < b.limit for b in bands) && levels_equal && !any(nonfinite .& active)
        pass &= ok
        finite_err = maximum(err[.!nonfinite])
        results[name] = Dict("max_error_active_ocean" => ocean_err, "max_error_away_from_singularities" => far_err,
                             "max_error_active_ocean_south_of_88N" => south_of_cap,
                             "max_error_active_ocean_north_of_88N" => cap,
                             "max_error_by_latitude_band" => Dict(b.name => b.error for b in bands),
                             "limit_by_latitude_band" => Dict(b.name => b.limit for b in bands),
                             "max_error_top_three_rows_active_ocean" => top_err, "max_error_south_row" => south_err,
                             "max_error_finite_cells" => finite_err, "nonfinite_cells" => count(nonfinite),
                             "nonfinite_cells_away_from_singularities" => count(nonfinite .& far),
                             "all_levels_equal" => levels_equal, "pass" => ok)
        println("SELF_TEST ", name, " ocean=", ocean_err, " top_rows_ocean=", top_err, " far=", far_err,
                " south_row=", south_err, " finite=", finite_err,
                " nonfinite=", count(nonfinite), " nonfinite_far=", count(nonfinite .& far), " pass=", ok)
    end
    return pass, results, lon, lat
end

function main(args)
    length(args) >= 6 || error("usage: NX NY NZ DEPTH_M OUT_DIR MASK_FILE ANNUAL_FILE...")
    nx, ny, nz = parse.(Int, args[1:3])
    depth = parse(Float64, args[4])
    out = args[5]
    ispath(out) && error("$out exists")
    files = args[7:end]
    grid = build_grid(nx, ny, nz, depth)
    active = NCDataset(ds -> ds["ocean_active_mask"][:, :] .== 1, args[6])
    u, v = velocity_fields(grid)
    println("U_NORTH_BC ", summary(u.boundary_conditions.north))
    println("V_NORTH_BC ", summary(v.boundary_conditions.north))
    pass, results, lon, lat = self_test(grid, active)
    pass || error("self-test failed: $results")
    mkpath(out)
    record = Dict{String,Any}("schema" => "readyesm_offline_geographic_velocity_v1", "observer_sha256" => OBSERVER_SHA256,
                              "observer_commit" => "019aae5", "script_sha256" => filesha(@__FILE__),
                              "grid" => Dict("size" => [nx, ny, nz], "depth_m" => depth, "halo" => [5, 5, 4]),
                              "self_test" => results, "outputs" => Dict{String,Any}())
    zc = [Float64(O.znode(1, 1, k, grid, O.Center(), O.Center(), O.Center())) for k in 1:nz]
    for path in files
        NCDataset(path) do ds
            # The grid must be the writer's grid: same centres in longitude, latitude and depth.
            dlon = maximum(abs.(mod.(ds["ocean_longitude"][:, :] .- lon .+ 180, 360) .- 180))
            dlat = maximum(abs.(ds["ocean_latitude"][:, :] .- lat))
            dz = maximum(abs.(ds["ocean_z"][:] .- zc))
            (dlon < 1e-3 && dlat < 1e-3 && dz < 1e-2) || error("$path: grid differs (lon $dlon, lat $dlat, z $dz)")
            uarr = Float32.(ds["cmip_native_u"][:, :, :])
            varr = Float32.(ds["cmip_native_v"][:, :, :])
            east, north = rotate(grid, uarr, varr)
            name = "geographic_velocity_" * replace(basename(path), "annual_" => "")
            target = joinpath(out, name)
            NCDataset(target, "c") do o
                defDim(o, "ocean_i", nx); defDim(o, "ocean_j", ny); defDim(o, "ocean_k", nz)
                for (var, values, std) in (("eastward_velocity", east, "eastward_sea_water_velocity"),
                                           ("northward_velocity", north, "northward_sea_water_velocity"))
                    v = defVar(o, var, values, ("ocean_i", "ocean_j", "ocean_k"))
                    v.attrib["units"] = "m s-1"
                    v.attrib["standard_name"] = std
                    v.attrib["grid_location"] = "horizontal_cell_center"
                    v.attrib["interpolation"] = _SURFACE_GEOGRAPHIC_VELOCITY_PROVENANCE
                end
                o.attrib["source_file"] = path
                o.attrib["source_sha256"] = filesha(path)
                o.attrib["observer_sha256"] = OBSERVER_SHA256
                o.attrib["script_sha256"] = filesha(@__FILE__)
            end
            record["outputs"][name] = Dict("source" => path, "source_sha256" => filesha(path),
                                           "sha256" => filesha(target), "grid_max_abs_diff" => [dlon, dlat, dz])
            println("ROTATED ", path, " -> ", target)
        end
    end
    open(joinpath(out, "record.toml"), "w") do io
        TOML.print(io, record)
    end
    println("ROTATE_NATIVE_VELOCITY_DONE")
end

main(ARGS)
