"""
Construct EVP auxiliaries on a latitude--longitude grid without evaluating the
constitutive kernels in geometrically invalid polar halo coordinates.

ClimaSeaIce's generic constructor evaluates viscosity and stress kernels from
`-halo+2` through `N+halo-1`. Compact latitude--longitude grids whose physical
domain approaches the poles extrapolate those halos beyond ±90 degrees, where
the signed face area becomes negative. CPU execution then throws when the EVP
relaxation parameter takes `sqrt(ζ cα Δt / (m A))`; CUDA silently creates NaNs
that are subsequently overwritten by the stress halo fill. Only interior
stresses are physically needed because `finalize_rheology!` fills their halos.
The tripolar constructor is deliberately untouched.
"""
function ClimaSeaIce.Rheologies.Auxiliaries(
    rheology::ClimaSeaIce.Rheologies.ElastoViscoPlasticRheology,
    grid::Oceananigans.LatitudeLongitudeGrid,
)
    R = ClimaSeaIce.Rheologies
    arch = Oceananigans.Architectures.architecture(grid)
    Nx, Ny, _ = size(grid)
    interior_parameters = Oceananigans.Utils.KernelParameters(1:Nx, 1:Ny)

    σ₁₁ = Oceananigans.Field{
        Oceananigans.Center,
        Oceananigans.Center,
        Nothing,
    }(grid)
    σ₂₂ = similar(σ₁₁)
    σ₁₂ = Oceananigans.Field{
        Oceananigans.Face,
        Oceananigans.Face,
        Nothing,
    }(grid)
    uⁿ = Oceananigans.Field{
        Oceananigans.Face,
        Oceananigans.Center,
        Nothing,
    }(grid)
    vⁿ = Oceananigans.Field{
        Oceananigans.Center,
        Oceananigans.Face,
        Nothing,
    }(grid)
    P = similar(σ₁₁)
    α = similar(σ₁₁)
    Δ = similar(σ₁₁)
    ζᶠᶠᶜ = similar(σ₁₂)
    ζᶜᶜᶜ = similar(σ₁₁)

    fill!(α, rheology.max_relaxation_parameter)

    viscosity_kernel = Oceananigans.Utils.configure_kernel(
        arch,
        grid,
        interior_parameters,
        R._compute_evp_viscosities!,
    )[1]
    stresses_kernel = Oceananigans.Utils.configure_kernel(
        arch,
        grid,
        interior_parameters,
        R._compute_evp_stresses!,
    )[1]

    initialization_parameters = Oceananigans.Utils.KernelParameters(
        size(P.data)[1:2],
        P.data.offsets[1:2],
    )
    initialization_kernel = Oceananigans.Utils.configure_kernel(
        arch,
        grid,
        initialization_parameters,
        R._initialize_evp_rhology!,
    )[1]

    fields = (; σ₁₁, σ₂₂, σ₁₂, ζᶠᶠᶜ, ζᶜᶜᶜ, Δ, α, uⁿ, vⁿ, P)
    kernels = (
        _viscosity_kernel! = viscosity_kernel,
        _stresses_kernel! = stresses_kernel,
        _initialize_rhology! = initialization_kernel,
    )
    return R.Auxiliaries(fields, kernels)
end

const CLIMASEAICE_LATLON_EVP_KERNEL_DOMAIN =
    "physical_interior_then_native_stress_halo_fill_v1"
