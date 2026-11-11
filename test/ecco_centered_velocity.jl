using Dates
using JLD2
using NCDatasets

# A lazy global source uses the actual 0.5-degree ECCO4 coordinate convention,
# without allocating or reading a 720x360x50 production volume in a unit test.
struct ECCOAffineSource{F} <: AbstractArray{Float32, 3}
    value::F
end
Base.size(::ECCOAffineSource) = (720, 360, 50)
Base.getindex(a::ECCOAffineSource, i::Int, j::Int, k::Int) =
    Float32(a.value(-179.75 + (i-1)/2, -89.75 + (j-1)/2))

function write_ecco_velocity_test_header(path, variable)
    NCDataset(path, "c") do ds
        for (name, n) in (("i", 720), ("j", 360), ("k", 50), ("time", 1))
            defDim(ds, name, n)
        end
        defVar(ds, "longitude", Float64, ("i",))[:] = collect(-179.75:0.5:179.75)
        defVar(ds, "latitude", Float64, ("j",))[:] = collect(-89.75:0.5:89.75)
        faces = ReadyESM.NumericalEarth.DataWrangling.z_interfaces(ReadyESM.NumericalEarth.ECCO4Monthly())
        defVar(ds, "Z", Float64, ("k",))[:] = reverse((faces[1:end-1] .+ faces[2:end]) ./ 2)
        # Chunked variable is left unallocated: this fixture tests the native
        # mapping kernel with the lazy source and cache loading separately.
        defVar(ds, variable, Float32, ("i", "j", "k", "time");
            chunksizes = (16, 16, 1, 1), deflatelevel = 1)
    end
end

@testset "ECCO4 centered velocity and cache contracts" begin
    O = ReadyESM.Oceananigans
    NE = ReadyESM.NumericalEarth
    DW = NE.DataWrangling
    region = NE.BoundingBox(longitude = (-160., -158.), latitude = (-30., -28.))
    eastward(λ, φ) = 2 + 0.01λ + 0.02φ
    northward(λ, φ) = -1 + 0.03λ - 0.01φ
    mktempdir() do dir
        sources = []
        for (name, variable, affine) in ((:u_velocity, "EVEL", eastward), (:v_velocity, "NVEL", northward))
            metadata = NE.Metadatum(name; dataset = NE.ECCO4Monthly(),
                date = DateTime(1993, 1, 1), dir, region)
            write_ecco_velocity_test_header(DW.metadata_path(metadata), variable)
            identity = ReadyESM._ecco4_velocity_cache_identity(metadata)
            field = O.CenterField(DW.native_grid(metadata, O.CPU()))
            ReadyESM._populate_ecco4_centered_velocity!(field, metadata, ECCOAffineSource(affine))
            for (i, j) in ((2, 2), (3, 4), (5, 5))
                λ, φ, z = O.Grids.node(i, j, 50, field.grid, O.Center(), O.Center(), O.Center())
                @test field[i, j, 50] ≈ affine(λ, φ) atol = 16eps(Float32)
            end
            @test_throws ArgumentError ReadyESM._populate_ecco4_centered_velocity!(
                O.XFaceField(field.grid), metadata, ECCOAffineSource(affine))
            legacy = DW.inpainted_metadata_path(metadata)
            JLD2.jldsave(legacy; data = fill(999f0, size(parent(field))), inpainting_maxiter = 5)
            @test_throws ErrorException ReadyESM._ecco4_centered_velocity_field(metadata)
            @test_throws ErrorException ReadyESM._read_ecco4_velocity_cache!(field, legacy, identity)
            cache = ReadyESM._ecco4_velocity_cache_path(metadata)
            @test cache != legacy
            ReadyESM._write_ecco4_velocity_cache(cache, field, identity, Dict("test" => true))
            loaded = ReadyESM._ecco4_centered_velocity_field(metadata)
            @test parent(loaded) == parent(field)
            @test_throws ErrorException ReadyESM._write_ecco4_velocity_cache(cache, field, identity, Dict())

            # Reject old schema, stale source hash, incompatible halo, corrupted
            # data and a changed source file even if array dimensions match.
            for key in ("schema", "source_sha256", "location", "halo")
                wrong = copy(identity)
                wrong[key] = "wrong"
                @test_throws ErrorException ReadyESM._read_ecco4_velocity_cache!(loaded, cache, wrong)
            end
            bad = joinpath(dir, variable * "_nan.jld2")
            data = copy(parent(field))
            data[1] = NaN32
            JLD2.jldsave(bad; identity, data)
            @test_throws ErrorException ReadyESM._read_ecco4_velocity_cache!(loaded, bad, identity)
            for data in (zeros(Float32, 1, 1, 1), Float64.(parent(field)), fill(1f6, size(parent(field))))
                rm(bad)
                JLD2.jldsave(bad; identity, data)
                @test_throws ErrorException ReadyESM._read_ecco4_velocity_cache!(loaded, bad, identity)
            end
            NCDataset(DW.metadata_path(metadata), "a") do ds
                ds.attrib["changed_source"] = "cache must be rebuilt"
            end
            @test_throws ErrorException ReadyESM._ecco4_centered_velocity_field(metadata)
            rm(cache)
            identity = ReadyESM._ecco4_velocity_cache_identity(metadata)
            ReadyESM._write_ecco4_velocity_cache(cache, field, identity, Dict("test" => true))
            push!(sources, ReadyESM.ECCO4CenteredVelocity(metadata))
        end

        grid = O.LatitudeLongitudeGrid(O.CPU(), Float64; size = (4, 4, 2),
            topology = (O.Bounded, O.Bounded, O.Bounded), halo = (3, 3, 3),
            longitude = (-159.75, -158.25), latitude = (-29.75, -28.25), z = (-15, -5))
        for (source, affine) in zip(sources, (eastward, northward))
            target = O.CenterField(grid)
            O.set!(target, source)
            for j in 1:4, i in 1:4
                λ, φ, z = O.Grids.node(i, j, 1, grid, O.Center(), O.Center(), O.Center())
                @test target[i, j, 1] ≈ affine(λ, φ) atol = 16eps(Float32)
            end
        end
        u, v = O.XFaceField(grid), O.YFaceField(grid)
        ReadyESM._set_extrinsic_ocean_velocity_fields!(u, v, grid, sources...)
        old = (u = O.XFaceField(grid), v = O.YFaceField(grid))
        O.Models.HydrostaticFreeSurfaceModels.set_from_extrinsic_velocities!(old, grid, sources...)
        for j in 2:3, i in 2:3
            λu, φu, _ = O.Grids.node(i, j, 1, grid, O.Face(), O.Center(), O.Center())
            λv, φv, _ = O.Grids.node(i, j, 1, grid, O.Center(), O.Face(), O.Center())
            @test u[i, j, 1] ≈ eastward(λu, φu) atol = 16eps(Float32)
            @test v[i, j, 1] ≈ northward(λv, φv) atol = 16eps(Float32)
            @test old.u[i, j, 1] - u[i, j, 1] ≈ 0.02 * 0.375 / 2 atol = 16eps(Float32)
            @test old.v[i, j, 1] - v[i, j, 1] ≈ 0.03 * 0.375 / 2 atol = 16eps(Float32)
        end
    end
end

@testset "Geographic vector rotation then face interpolation" begin
    O = ReadyESM.Oceananigans
    grid = O.TripolarGrid(O.CPU(), Float64; size = (16, 12, 2), halo = (3, 3, 3), z = (-100, 0))
    u, v = O.XFaceField(grid), O.YFaceField(grid)
    ReadyESM._set_extrinsic_ocean_velocity_fields!(u, v, grid, 2.0, 0.5)
    ur(i, j) = 2cos(O.Operators.rotation_angle(i, j, grid)) - 0.5sin(O.Operators.rotation_angle(i, j, grid))
    vr(i, j) = 2sin(O.Operators.rotation_angle(i, j, grid)) + 0.5cos(O.Operators.rotation_angle(i, j, grid))
    for (i, j) in ((3, 3), (7, 8), (5, 10), (12, 10))
        @test u[i, j, 1] ≈ (ur(i-1, j) + ur(i, j))/2 atol = 1e-12
        @test v[i, j, 1] ≈ (vr(i, j-1) + vr(i, j))/2 atol = 1e-12
        @test ur(i, j)^2 + vr(i, j)^2 ≈ 4.25 atol = 1e-12
    end
    # The last centered row is a duplicated fold: its second half represents
    # the same points in the opposite native frame. Derive the reference from
    # the signed geometric mapping, not from a similarly halo-filled field.
    folded(f, i, j) = j == 12 && i > 8 ? -f(17-i, j) : f(i, j)
    for i in 1:16
        @test u[i, 12, 1] ≈ (folded(ur, mod1(i-1, 16), 12) + folded(ur, i, 12))/2 atol = 1e-12
        @test v[i, 12, 1] ≈ (vr(i, 11) + folded(vr, i, 12))/2 atol = 1e-12
    end
    @test abs(u[1, 12, 1]) < 1e-12
    @test abs(u[9, 12, 1]) < 1e-12
end

@testset "Bounded ECCO velocity inpainting" begin
    O = ReadyESM.Oceananigans
    grid = O.RectilinearGrid(O.CPU(), Float64; size = (5, 5, 2),
        topology = (O.Periodic, O.Bounded, O.Bounded), extent = (5, 5, 2))
    for donor_case in (:zero, :canceling, :nonzero)
        field = O.CenterField(grid)
        mask = O.CenterField(grid, Bool)
        O.set!(field, donor_case == :nonzero ? 7 : 0)
        mask[3, 3, 2] = true
        if donor_case == :canceling
            field[2, 3, 2], field[4, 3, 2] = -2, 2
        end
        statistics = ReadyESM._inpaint_ecco4_velocity!(field, mask)
        @test field[3, 3, 2] == (donor_case == :nonzero ? 7 : 0)
        @test statistics["horizontal_filled_cells"] == 1
        @test statistics["zero_fallback_cells"] == 0
        @test statistics["zero_mean_filled_cells"] == (donor_case == :nonzero ? 0 : 1)
        @test statistics["source_valid_max_absolute_change"] == 0
    end
    field, mask = O.CenterField(grid), O.CenterField(grid, Bool)
    O.set!(field, 7)
    O.interior(mask)[:, :, 1] .= true
    statistics = ReadyESM._inpaint_ecco4_velocity!(field, mask)
    @test all(==(7), O.interior(field))
    @test statistics["downward_filled_cells"] == 25
    @test statistics["zero_fallback_cells"] == 0
    O.set!(mask, true)
    statistics = ReadyESM._inpaint_ecco4_velocity!(field, mask)
    @test all(iszero, O.interior(field))
    @test statistics["horizontal_passes"] == 0
    @test statistics["zero_fallback_cells"] == 50
end
