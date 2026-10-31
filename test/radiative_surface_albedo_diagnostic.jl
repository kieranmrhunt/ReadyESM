module RadiativeSurfaceAlbedoDiagnosticTests
using ReadyESM, Test
const R = ReadyESM

@testset "Radiative surface albedo diagnostic" begin
    # Band-invariant surface (default): band 1 exactly, and the old refusals stay.
    invariant = repeat([0.1 0.4 0.85], 14, 1)
    @test R._radiative_surface_albedo_from_bands(invariant, copy(invariant), 3) == [0.1, 0.4, 0.85]
    @test R._radiative_surface_albedo_from_bands(invariant, copy(invariant), 3; spectral = true) ≈ [0.1, 0.4, 0.85]
    direct = [0.05 0.40; 0.07 0.40; 0.30 0.40; 0.34 0.40]
    diffuse = [0.06 0.40; 0.08 0.40; 0.32 0.40; 0.36 0.40]
    @test_throws ErrorException R._radiative_surface_albedo_from_bands(direct, diffuse, 2)
    banded = [0.1 0.4; 0.3 0.4]
    @test_throws ErrorException R._radiative_surface_albedo_from_bands(banded, copy(banded), 2)
    # Spectral surface: unweighted mean over bands and over direct and diffuse.
    @test R._radiative_surface_albedo_from_bands(direct, diffuse, 2; spectral = true) ≈ [0.1975, 0.40]
    # Shape, range and finiteness checks hold in both modes.
    @test_throws DimensionMismatch R._radiative_surface_albedo_from_bands(direct, diffuse, 3; spectral = true)
    @test_throws DimensionMismatch R._radiative_surface_albedo_from_bands(direct, diffuse[:, 1:1], 2; spectral = true)
    @test_throws ErrorException R._radiative_surface_albedo_from_bands(direct .+ 1, diffuse .+ 1, 2; spectral = true)
    bad = copy(direct); bad[1, 1] = NaN
    @test_throws ErrorException R._radiative_surface_albedo_from_bands(bad, diffuse, 2; spectral = true)
    # Out-of-range band values must be refused even when their mean is in range (Codex review, 08:22).
    cancelling = [-0.2 0.4; 1.2 0.4]
    @test_throws ErrorException R._radiative_surface_albedo_from_bands(cancelling, copy(cancelling), 2; spectral = true)
end
end
