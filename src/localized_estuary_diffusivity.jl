"""
    LocalizedEstuaryDiffusivity(vertical_diffusivity, horizontal_diffusivity,
                                tracer_selector=AllLocalizedEstuaryTracers())

An anisotropic tracer-only Oceananigans closure for unresolved river mouths.
`vertical_diffusivity` is applied with vertically implicit time discretization;
`horizontal_diffusivity` is applied by Oceananigans' ordinary conservative
flux-form horizontal operator. Both coefficients are cell-centred fields.
`tracer_selector` makes the versioned river-plume repair either reproduce the
legacy all-tracer behavior or act only on the salinity tracer.

Keeping the two directional coefficients in one closure is important for the
pinned Oceananigans version: the production ocean already has four closures,
and its tuple kernels are explicitly unrolled only through five entries. Two
separate river-mouth closures therefore enter an allocating six-entry fallback
on GPU. This type preserves the same two diffusion operators while keeping the
production tuple at five entries.
"""
abstract type AbstractLocalizedEstuaryTracerSelector end

struct AllLocalizedEstuaryTracers <:
       AbstractLocalizedEstuaryTracerSelector end

struct LocalizedEstuaryTracerIndex{I} <:
       AbstractLocalizedEstuaryTracerSelector end

@inline _localized_estuary_tracer_active(
    ::AllLocalizedEstuaryTracers,
    ::Val,
) = true

@inline _localized_estuary_tracer_active(
    ::LocalizedEstuaryTracerIndex{I},
    ::Val{I},
) where {I} = true

@inline _localized_estuary_tracer_active(
    ::LocalizedEstuaryTracerIndex,
    ::Val,
) = false

struct LocalizedEstuaryDiffusivity{V, H, S} <:
       Oceananigans.TurbulenceClosures.AbstractScalarDiffusivity{
           Oceananigans.TimeSteppers.VerticallyImplicitTimeDiscretization,
           Oceananigans.TurbulenceClosures.ThreeDimensionalFormulation,
           1,
    }
    vertical_diffusivity::V
    horizontal_diffusivity::H
    tracer_selector::S
end

LocalizedEstuaryDiffusivity(vertical_diffusivity, horizontal_diffusivity) =
    LocalizedEstuaryDiffusivity(
        vertical_diffusivity,
        horizontal_diffusivity,
        AllLocalizedEstuaryTracers(),
    )

function Base.summary(::LocalizedEstuaryDiffusivity)
    return "LocalizedEstuaryDiffusivity{VerticallyImplicitVertical," *
           "ExplicitHorizontal}"
end

Base.show(io::IO, closure::LocalizedEstuaryDiffusivity) =
    print(io, summary(closure))

function Adapt.adapt_structure(to, closure::LocalizedEstuaryDiffusivity)
    return LocalizedEstuaryDiffusivity(
        Adapt.adapt(to, closure.vertical_diffusivity),
        Adapt.adapt(to, closure.horizontal_diffusivity),
        closure.tracer_selector,
    )
end

function Oceananigans.Architectures.on_architecture(
    architecture,
    closure::LocalizedEstuaryDiffusivity,
)
    return LocalizedEstuaryDiffusivity(
        Oceananigans.on_architecture(
            architecture,
            closure.vertical_diffusivity,
        ),
        Oceananigans.on_architecture(
            architecture,
            closure.horizontal_diffusivity,
        ),
        closure.tracer_selector,
    )
end

# This closure transports tracers only. The ThreeDimensionalFormulation base
# methods also query viscosity while forming momentum tendencies; an exact zero
# retains their no-momentum effect without introducing special momentum kernels.
@inline function Oceananigans.TurbulenceClosures.viscosity(
    closure::LocalizedEstuaryDiffusivity,
    closure_fields,
)
    return zero(eltype(closure.vertical_diffusivity))
end

# Match ScalarDiffusivity's interpolation of cell-centred coefficient fields to
# the three tracer-flux faces, but select the intended directional field.
@inline function Oceananigans.TurbulenceClosures.κhᶠᶜᶜ(
    i,
    j,
    k,
    grid,
    closure::LocalizedEstuaryDiffusivity,
    closure_fields,
    tracer_id,
    clock,
    fields,
)
    if _localized_estuary_tracer_active(
        closure.tracer_selector,
        tracer_id,
    )
        return Oceananigans.Operators.ℑxᶠᵃᵃ(
            i,
            j,
            k,
            grid,
            closure.horizontal_diffusivity,
        )
    else
        return zero(eltype(closure.horizontal_diffusivity))
    end
end

@inline function Oceananigans.TurbulenceClosures.κhᶜᶠᶜ(
    i,
    j,
    k,
    grid,
    closure::LocalizedEstuaryDiffusivity,
    closure_fields,
    tracer_id,
    clock,
    fields,
)
    if _localized_estuary_tracer_active(
        closure.tracer_selector,
        tracer_id,
    )
        return Oceananigans.Operators.ℑyᵃᶠᵃ(
            i,
            j,
            k,
            grid,
            closure.horizontal_diffusivity,
        )
    else
        return zero(eltype(closure.horizontal_diffusivity))
    end
end

@inline function Oceananigans.TurbulenceClosures.κzᶜᶜᶠ(
    i,
    j,
    k,
    grid,
    closure::LocalizedEstuaryDiffusivity,
    closure_fields,
    tracer_id,
    clock,
    fields,
)
    if _localized_estuary_tracer_active(
        closure.tracer_selector,
        tracer_id,
    )
        return Oceananigans.Operators.ℑzᵃᵃᶠ(
            i,
            j,
            k,
            grid,
            closure.vertical_diffusivity,
        )
    else
        return zero(eltype(closure.vertical_diffusivity))
    end
end
