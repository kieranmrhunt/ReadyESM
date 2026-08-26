const _NUMERICAL_EARTH_SPEEDY_EXTENSION =
    Base.get_extension(NumericalEarth, :NumericalEarthSpeedyWeatherExt)

isnothing(_NUMERICAL_EARTH_SPEEDY_EXTENSION) && error(
    "NumericalEarthSpeedyWeatherExt was not loaded after importing NumericalEarth and SpeedyWeather",
)

const _ReadyESMSpeedySimulation =
    _NUMERICAL_EARTH_SPEEDY_EXTENSION.SpeedySimulation

"""
Rotate geographic velocity components on a tripolar grid using vector-aware
centered-field fold halos.

Oceananigans' generic initializer creates unnamed centered scratch fields.
Those fields receive the scalar zipper condition (`+1`), even though their
contents become intrinsic vector components before their north halos are
filled.  The staggered interpolation then reads a reflected rather than
sign-reversed component at the tripolar fold.  NumericalEarth already exposes
the boundary-condition constructor used by its other tripolar vector fields;
use it here while otherwise preserving Oceananigans' reference algorithm.
"""
function Oceananigans.Models.HydrostaticFreeSurfaceModels.set_from_extrinsic_velocities!(
    velocities,
    grid::Oceananigans.OrthogonalSphericalShellGrids.TripolarGridOfSomeKind,
    u,
    v,
)
    arch = Oceananigans.architecture(grid)
    vector_bcs = NumericalEarth.EarthSystemModels.InterfaceComputations.vector_component_boundary_conditions(
        grid,
        (Oceananigans.Center(), Oceananigans.Center(), Oceananigans.Center()),
    )
    uᶜᶜᶜ = Oceananigans.CenterField(grid; boundary_conditions = vector_bcs)
    vᶜᶜᶜ = Oceananigans.CenterField(grid; boundary_conditions = vector_bcs)
    u isa Oceananigans.Fields.ZeroField || Oceananigans.set!(uᶜᶜᶜ, u)
    v isa Oceananigans.Fields.ZeroField || Oceananigans.set!(vᶜᶜᶜ, v)
    Oceananigans.Utils.launch!(
        arch,
        grid,
        :xyz,
        Oceananigans.Models.HydrostaticFreeSurfaceModels._rotate_velocities!,
        uᶜᶜᶜ,
        vᶜᶜᶜ,
        grid,
    )
    Oceananigans.fill_halo_regions!(uᶜᶜᶜ)
    Oceananigans.fill_halo_regions!(vᶜᶜᶜ)
    Oceananigans.Utils.launch!(
        arch,
        grid,
        :xyz,
        Oceananigans.Models.HydrostaticFreeSurfaceModels._interpolate_velocities!,
        velocities.u,
        velocities.v,
        grid,
        uᶜᶜᶜ,
        vᶜᶜᶜ,
    )
    return nothing
end

_profile_regridding_construction() =
    lowercase(get(ENV, "READYESM_PROFILE_REGRIDDING_CONSTRUCTION", "false")) in
    ("1", "true", "yes", "on")

function _profile_regridding_phase(function_, architecture, label)
    _profile_regridding_construction() || return function_()
    Oceananigans.Architectures.synchronize(architecture)
    start = time_ns()
    value = function_()
    Oceananigans.Architectures.synchronize(architecture)
    seconds = (time_ns() - start) / 1e9
    println("REGRIDDING_CONSTRUCTION_PHASE label=$label seconds=$seconds")
    flush(stdout)
    return value
end

"""
Return a CPU, Float64 copy of the horizontal Oceananigans grid for robust
spherical-polygon intersection. Immersed bathymetry does not change the
horizontal cell geometry, so only its underlying grid is needed here.
"""
function _float64_regridding_geometry_grid(exchange_grid)
    underlying_grid = Oceananigans.ImmersedBoundaries.underlying_grid(exchange_grid)
    cpu_grid = Oceananigans.on_architecture(Oceananigans.CPU(), underlying_grid)

    # Oceananigans does not currently define `with_number_type` for an
    # OrthogonalSphericalShellGrid. Rebuild a geometrically identical,
    # one-layer tripolar grid explicitly. Only the horizontal coordinates are
    # consumed by ConservativeRegridding, so neither the live vertical grid nor
    # the immersed bathymetry belongs in this temporary geometry object.
    mapping = hasproperty(cpu_grid, :conformal_mapping) ?
              getproperty(cpu_grid, :conformal_mapping) : nothing
    if mapping isa Oceananigans.OrthogonalSphericalShellGrids.Tripolar
        Nx, Ny, _ = size(cpu_grid)
        Hx, Hy, _ = Oceananigans.halo_size(cpu_grid)
        fold_topology = Oceananigans.topology(cpu_grid)[2]

        return Oceananigans.TripolarGrid(
            Oceananigans.CPU(),
            Float64;
            size = (Nx, Ny, 1),
            halo = (Hx, Hy, 1),
            z = (0.0, 1.0),
            radius = Float64(cpu_grid.radius),
            southernmost_latitude = Float64(mapping.southernmost_latitude),
            north_poles_latitude = Float64(mapping.north_poles_latitude),
            first_pole_longitude = Float64(mapping.first_pole_longitude),
            fold_topology,
        )
    end

    return Oceananigans.Grids.with_number_type(Float64, cpu_grid)
end

"""
Convert a CPU regridder, including its sparse matrix and work vectors, to
Float32 before moving it to CUDA. ConservativeRegridding defaults to Float64;
leaving those arrays unchanged makes cuSPARSE reject Float32 model fields.
"""
function _float32_regridder(regridder)
    return ConservativeRegridding.Regridder(
        Float32.(regridder.intersections),
        Float32.(regridder.dst_areas),
        Float32.(regridder.src_areas),
        Float32.(regridder.dst_temp),
        Float32.(regridder.src_temp),
    )
end

"""
Materialize one GPU regridding direction as an explicit CSR matrix.

Oceananigans' generic architecture transfer converts a CPU `SparseMatrixCSC`
to `CuSparseMatrixCSC`. Multiplication by that column-oriented matrix scatters
contributions into destination rows and was observed to produce different
Float32 bits on repeated identical full-grid applications. CSR instead assigns
one output row to each sparse dot product. The reverse direction must be
materialized independently as CSR(Aᵀ); a lazy transpose would reintroduce a
transposed sparse multiplication path and shared work buffers.
"""
function _explicit_row_regridder(architecture, regridder)
    architecture isa Oceananigans.CPU && return regridder
    to_architecture(values) = Oceananigans.on_architecture(
        architecture,
        values,
    )
    intersections = CUDA.CUSPARSE.CuSparseMatrixCSR{Float32}(
        regridder.intersections,
    )
    return ConservativeRegridding.Regridder(
        intersections,
        to_architecture(regridder.dst_areas),
        to_architecture(regridder.src_areas),
        to_architecture(regridder.dst_temp),
        to_architecture(regridder.src_temp),
    )
end

"""
Apply a non-transposed CUDA CSR operator with cuSPARSE's bitwise-deterministic
SpMV algorithm.

`LinearAlgebra.mul!` selects `CUSPARSE_SPMV_ALG_DEFAULT`, which maps CSR to
`CUSPARSE_SPMV_CSR_ALG1`; NVIDIA documents that algorithm as potentially
different between identical runs.  Algorithm 2 is deterministic for a
non-transposed CSR operator.  Clear the destination first and call SpMV with
`beta = 1`: this is algebraically `A * source + 0`, and avoids cuSPARSE's
documented Compute Sanitizer false-race optimization for `beta = 0`.
"""
function _deterministic_csr_spmv!(
    destination::CUDA.CuVector{T},
    intersections::CUDA.CUSPARSE.CuSparseMatrixCSR,
    source::CUDA.CuVector{T},
) where {T}
    fill!(destination, zero(T))
    CUDA.CUSPARSE.mv!(
        'N',
        one(T),
        intersections,
        source,
        one(T),
        destination,
        'O',
        CUDA.CUSPARSE.CUSPARSE_SPMV_CSR_ALG2,
    )
    return destination
end

function ConservativeRegridding.perform_regridding!(
    destination::CUDA.CuVector{T},
    regridder::ConservativeRegridding.Regridder{W},
    source::CUDA.CuVector{T};
    kwargs...,
) where {T, W <: CUDA.CUSPARSE.CuSparseMatrixCSR}
    _deterministic_csr_spmv!(destination, regridder.intersections, source)
    return destination
end

"""
Float32 Oceananigans specialization of NumericalEarth's SpeedyWeather state
exchanger. Compute cell intersections using Float64 coordinates on the CPU,
then use Float32 weights and buffers for the live CPU/GPU model exchange.

This avoids a mixed Float64/Float32 GeometryOps clipping failure during model
construction and preserves type-matched Float32 sparse matvecs on CUDA.
"""
function NumericalEarth.EarthSystemModels.InterfaceComputations.ComponentExchanger(
    atmosphere::_ReadyESMSpeedySimulation,
    exchange_grid::Oceananigans.Grids.AbstractGrid{Float32};
    correction = nothing,
)
    arch = Oceananigans.architecture(exchange_grid)
    atmosphere_architecture = atmosphere.model.spectral_grid.architecture
    if (arch isa Oceananigans.CPU) !=
       (atmosphere_architecture isa SpeedyWeather.CPU)
        error(
            "The exchange grid is on $arch but the SpeedyWeather atmosphere is on " *
            "$atmosphere_architecture. Both components must use the same architecture.",
        )
    end

    spectral_grid = atmosphere.model.spectral_grid.grid
    cpu_spectral_grid = _profile_regridding_phase(
        arch,
        "cpu_spectral_grid",
    ) do
        SpeedyWeather.on_architecture(SpeedyWeather.CPU(), spectral_grid)
    end
    geometry_grid = _profile_regridding_phase(
        arch,
        "float64_ocean_geometry",
    ) do
        _float64_regridding_geometry_grid(exchange_grid)
    end
    manifold = _profile_regridding_phase(arch, "best_manifold") do
        ConservativeRegridding.GOCore.best_manifold(geometry_grid)
    end

    cpu_from_atmosphere = _profile_regridding_phase(
        arch,
        "spherical_polygon_intersections",
    ) do
        ConservativeRegridding.Regridder(
            manifold,
            geometry_grid,
            cpu_spectral_grid,
        )
    end
    cpu_from_atmosphere = _profile_regridding_phase(
        arch,
        "float32_cpu_regridder",
    ) do
        _float32_regridder(cpu_from_atmosphere)
    end

    # Reuse the exact reciprocal CPU geometry, but on GPU materialize A and Aᵀ
    # independently in row-oriented storage. This avoids both nondeterministic
    # CSC scatter accumulation and a lazy transposed sparse-matvec path. CPU
    # execution retains ConservativeRegridding's ordinary reciprocal view.
    from_atmosphere = _profile_regridding_phase(
        arch,
        "forward_row_regridder",
    ) do
        _explicit_row_regridder(arch, cpu_from_atmosphere)
    end
    to_atmosphere = _profile_regridding_phase(
        arch,
        "reverse_row_regridder",
    ) do
        _explicit_row_regridder(arch, transpose(cpu_from_atmosphere))
    end
    regridder = (; to_atmosphere, from_atmosphere)

    Field2D = Oceananigans.Field{
        Oceananigans.Center,
        Oceananigans.Center,
        Nothing,
    }
    velocity_bcs = NumericalEarth.EarthSystemModels.InterfaceComputations.vector_component_boundary_conditions(
        exchange_grid,
        (Oceananigans.Center(), Oceananigans.Center(), nothing),
    )
    state = _profile_regridding_phase(arch, "exchange_state_fields") do
        (
            u = Field2D(exchange_grid; boundary_conditions = velocity_bcs),
            v = Field2D(exchange_grid; boundary_conditions = velocity_bcs),
            T = Field2D(exchange_grid),
            p = Field2D(exchange_grid),
            q = Field2D(exchange_grid),
            ℐꜜˢʷ = Field2D(exchange_grid),
            ℐꜜˡʷ = Field2D(exchange_grid),
            Jʳⁿ = Field2D(exchange_grid),
            Jˢⁿ = Field2D(exchange_grid),
            tmp = Field2D(exchange_grid),
        )
    end

    correction = _profile_regridding_phase(arch, "materialize_correction") do
        NumericalEarth.EarthSystemModels.InterfaceComputations.materialize_correction(
            correction,
            exchange_grid,
            atmosphere,
        )
    end
    return NumericalEarth.EarthSystemModels.InterfaceComputations.ComponentExchanger(
        state,
        regridder,
        correction,
    )
end

const _FRESHWATER_DENSITY_KG_M3 = 1000f0
# Key the cache by SpeedyWeather's mutable clock rather than by the immutable
# atmosphere wrapper. Weak ownership lets a completed coupled model release its
# exchange fields in interactive or matched multi-case processes.
const _INTENSIVE_REGRID_CACHE = WeakKeyDict{Any, Any}()

KernelAbstractions.@kernel function _mirror_right_center_fold_partners_kernel!(
    field,
    Nx,
    Ny,
    Nh,
    Nq,
)
    index = @index(Global)
    primary = ifelse(index <= Nq, index, index + Nh - Nq)
    partner = Nx + 1 - primary
    @inbounds field[partner, Ny, 1] = field[primary, Ny, 1]
end

"""
Fill the duplicate partner slots in a RightCenterFolded tripolar fold row.

ConservativeRegridding assigns zero area to these duplicate storage slots and
expects its Oceananigans field finalizer to mirror the unique physical cells
after normalization. NumericalEarth bypasses that finalizer by regridding into
`vec(interior(field))`, leaving `0 / 0 = NaN`. Mirroring must happen before
fold halo filling, and before geographic winds are rotated into grid axes, or
the invalid partner values can contaminate their physical counterparts.
"""
function _mirror_right_center_fold_partners!(field)
    grid = Oceananigans.ImmersedBoundaries.underlying_grid(field.grid)
    Oceananigans.topology(grid)[2] == Oceananigans.RightCenterFolded ||
        return field
    Nx, Ny, _ = size(grid)
    Nx % 4 == 0 || error(
        "RightCenterFolded conservative exchange requires longitude size divisible by 4",
    )
    Nh = Nx ÷ 2
    Nq = Nx ÷ 4
    interior = Oceananigans.interior(field)
    backend = KernelAbstractions.get_backend(interior)
    _mirror_right_center_fold_partners_kernel!(backend)(
        interior,
        Nx,
        Ny,
        Nh,
        Nq;
        ndrange = 2Nq,
    )
    return field
end

KernelAbstractions.@kernel function _wet_surface_mask_kernel!(wet_mask, grid)
    i, j = @index(Global, NTuple)
    k = size(grid, 3)
    inactive = Oceananigans.Grids.inactive_node(
        i,
        j,
        k,
        grid,
        Oceananigans.Center(),
        Oceananigans.Center(),
        Oceananigans.Center(),
    )
    @inbounds wet_mask[i, j, 1] = ifelse(inactive, zero(grid), one(grid))
end

KernelAbstractions.@kernel function _masked_mixed_surface_temperature_kernel!(
    destination,
    wet_mask,
    ocean_temperature,
    ice_temperature,
    ice_concentration,
    kelvin_offset,
)
    i, j = @index(Global, NTuple)
    @inbounds begin
        wet = wet_mask[i, j, 1] > zero(kelvin_offset)
        ice = ice_concentration[i, j, 1]
        mixed_temperature = ocean_temperature[i, j, 1] * (one(ice) - ice) +
                            ice * ice_temperature[i, j, 1] + kelvin_offset
        destination[i, j, 1] = ifelse(
            wet,
            mixed_temperature,
            zero(mixed_temperature),
        )
    end
end

KernelAbstractions.@kernel function _masked_intensive_source_kernel!(
    destination,
    wet_mask,
    source,
)
    i, j = @index(Global, NTuple)
    @inbounds begin
        value = source[i, j, 1]
        destination[i, j, 1] = ifelse(
            wet_mask[i, j, 1] > zero(value),
            value,
            zero(value),
        )
    end
end

@inline function _mixed_surface_flux(
    wet,
    ocean_flux,
    ice_flux,
    ice_concentration,
)
    if !wet
        return zero(ocean_flux)
    elseif ice_concentration <= zero(ice_concentration)
        # Do not form `0 * ice_flux`: the ice-side diagnostic is allowed to be
        # undefined where there is no sea ice.
        return ocean_flux
    elseif ice_concentration >= one(ice_concentration)
        # Likewise, do not evaluate an inactive ocean-side flux under complete
        # ice cover.
        return ice_flux
    end
    return ocean_flux * (one(ice_concentration) - ice_concentration) +
           ice_concentration * ice_flux
end

KernelAbstractions.@kernel function _masked_mixed_surface_flux_kernel!(
    destination,
    wet_mask,
    ocean_flux,
    ice_flux,
    ice_concentration,
)
    i, j = @index(Global, NTuple)
    @inbounds begin
        wet = wet_mask[i, j, 1] > zero(wet_mask[i, j, 1])
        destination[i, j, 1] = _mixed_surface_flux(
            wet,
            ocean_flux[i, j, 1],
            ice_flux[i, j, 1],
            ice_concentration[i, j, 1],
        )
    end
end

@inline function _coverage_normalized_intensive(
    numerator,
    coverage,
    fallback,
    minimum_coverage,
)
    return ifelse(
        coverage > minimum_coverage,
        numerator / max(coverage, minimum_coverage),
        fallback,
    )
end

function _intensive_regrid_cache!(coupled_model, atmosphere, tmp, to_atmosphere)
    cache_owner = atmosphere.variables.prognostic.clock
    return get!(_INTENSIVE_REGRID_CACHE, cache_owner) do
        grid = coupled_model.interfaces.exchanger.grid
        wet_mask = Oceananigans.Field{
            Oceananigans.Center,
            Oceananigans.Center,
            Nothing,
        }(grid)
        wet_interior = Oceananigans.interior(wet_mask)
        backend = KernelAbstractions.get_backend(wet_interior)
        _wet_surface_mask_kernel!(backend)(
            wet_interior,
            grid;
            ndrange = (size(grid, 1), size(grid, 2)),
        )
        KernelAbstractions.synchronize(backend)

        sst = atmosphere.variables.prognostic.ocean.sea_surface_temperature.data
        coverage = similar(sst)
        numerator = similar(sst)
        ConservativeRegridding.regrid!(coverage, to_atmosphere, wet_mask)
        KernelAbstractions.synchronize(backend)
        all(isfinite, Array(coverage)) || error(
            "ocean-to-atmosphere wet-area coverage contains non-finite values",
        )
        minimum(coverage) >= zero(eltype(coverage)) || error(
            "ocean-to-atmosphere wet-area coverage contains negative values",
        )
        maximum(coverage) <= one(eltype(coverage)) + 32eps(eltype(coverage)) ||
            error("ocean-to-atmosphere wet-area coverage exceeds unity")
        return (; wet_mask, coverage, numerator, backend)
    end
end

function _regrid_surface_intensives!(coupled_model, atmosphere)
    exchanger = coupled_model.interfaces.exchanger.atmosphere
    to_atmosphere = exchanger.regridder.to_atmosphere
    tmp = exchanger.state.tmp
    cache = _intensive_regrid_cache!(
        coupled_model,
        atmosphere,
        tmp,
        to_atmosphere,
    )
    grid = coupled_model.interfaces.exchanger.grid
    ocean_temperature =
        coupled_model.interfaces.atmosphere_ocean_interface.temperature
    ice_temperature =
        coupled_model.interfaces.atmosphere_sea_ice_interface.temperature
    ice_concentration =
        NumericalEarth.EarthSystemModels.sea_ice_concentration(coupled_model.sea_ice)
    sst = atmosphere.variables.prognostic.ocean.sea_surface_temperature.data
    ice_on_atmosphere_grid =
        atmosphere.variables.prognostic.ocean.sea_ice_concentration.data
    minimum_coverage = convert(eltype(sst), 1e-6)

    _masked_mixed_surface_temperature_kernel!(cache.backend)(
        Oceananigans.interior(tmp),
        Oceananigans.interior(cache.wet_mask),
        ocean_temperature,
        ice_temperature,
        ice_concentration,
        convert(eltype(sst), 273.15);
        ndrange = (size(grid, 1), size(grid, 2)),
    )
    KernelAbstractions.synchronize(cache.backend)
    ConservativeRegridding.regrid!(cache.numerator, to_atmosphere, tmp)
    sst .= _coverage_normalized_intensive.(
        cache.numerator,
        cache.coverage,
        convert(eltype(sst), 273.15),
        minimum_coverage,
    )

    _masked_intensive_source_kernel!(cache.backend)(
        Oceananigans.interior(tmp),
        Oceananigans.interior(cache.wet_mask),
        ice_concentration;
        ndrange = (size(grid, 1), size(grid, 2)),
    )
    KernelAbstractions.synchronize(cache.backend)
    ConservativeRegridding.regrid!(cache.numerator, to_atmosphere, tmp)
    ice_on_atmosphere_grid .= clamp.(
        _coverage_normalized_intensive.(
            cache.numerator,
            cache.coverage,
            zero(eltype(ice_on_atmosphere_grid)),
            minimum_coverage,
        ),
        zero(eltype(ice_on_atmosphere_grid)),
        one(eltype(ice_on_atmosphere_grid)),
    )
    return nothing
end

"""
Conservatively map turbulent ocean/sea-ice fluxes after masking inactive ocean
cells and selecting the physically active component before doing arithmetic.

NumericalEarth's generic SpeedyWeather exchange forms `(1-A) * ocean + A * ice`
on every exchange cell. Inactive immersed cells and absent surface components
may legitimately carry undefined diagnostics, so eager arithmetic such as
`0 * NaN` contaminates the sparse regrid. These are extensive whole-cell
fluxes: dry source cells are zeroed, but represented wet-area coverage is not
normalized as it is for SST and ice concentration.
"""
function _regrid_surface_fluxes!(coupled_model, atmosphere)
    exchanger = coupled_model.interfaces.exchanger.atmosphere
    to_atmosphere = exchanger.regridder.to_atmosphere
    tmp = exchanger.state.tmp
    cache = _intensive_regrid_cache!(
        coupled_model,
        atmosphere,
        tmp,
        to_atmosphere,
    )
    grid = coupled_model.interfaces.exchanger.grid
    atmosphere_ocean =
        coupled_model.interfaces.atmosphere_ocean_interface.fluxes
    atmosphere_ice =
        coupled_model.interfaces.atmosphere_sea_ice_interface.fluxes
    ice_concentration =
        NumericalEarth.EarthSystemModels.sea_ice_concentration(coupled_model.sea_ice)
    diagnosed_ocean = atmosphere.variables.parameterizations.ocean

    for (destination, ocean_flux, ice_flux) in (
        (
            diagnosed_ocean.sensible_heat_flux.data,
            atmosphere_ocean.sensible_heat,
            atmosphere_ice.sensible_heat,
        ),
        (
            diagnosed_ocean.surface_humidity_flux.data,
            atmosphere_ocean.water_vapor,
            atmosphere_ice.water_vapor,
        ),
    )
        _masked_mixed_surface_flux_kernel!(cache.backend)(
            Oceananigans.interior(tmp),
            Oceananigans.interior(cache.wet_mask),
            ocean_flux,
            ice_flux,
            ice_concentration;
            ndrange = (size(grid, 1), size(grid, 2)),
        )
        KernelAbstractions.synchronize(cache.backend)
        ConservativeRegridding.regrid!(destination, to_atmosphere, tmp)
    end
    return nothing
end

"""
Convert SpeedyWeather precipitation depth rates (m s⁻¹) to the freshwater
mass-flux convention used by NumericalEarth (kg m⁻² s⁻¹).

NumericalEarth divides these fields by the ocean reference density when it
assembles the free-surface flux.  Passing SpeedyWeather's depth rate through
unchanged therefore suppresses rain and snow freshwater input by a factor of
approximately 1000.
"""
function _precipitation_depth_to_mass_flux!(rain, snow)
    # Oceananigans Field data are OffsetArrays. Broadcasting directly into a
    # GPU-backed OffsetArray can select Base's scalar host loop instead of the
    # CuArray broadcast backend. Operate on the parent storage so CPU arrays
    # retain identical behavior and CUDA launches a device kernel. Include the
    # halo storage because NumericalEarth has already filled it in depth-rate
    # units before this conversion.
    parent(rain) .*= _FRESHWATER_DENSITY_KG_M3
    parent(snow) .*= _FRESHWATER_DENSITY_KG_M3
    return nothing
end

function NumericalEarth.EarthSystemModels.interpolate_state!(
    exchanger,
    exchange_grid::Oceananigans.Grids.AbstractGrid{Float32},
    atmosphere::_ReadyESMSpeedySimulation,
    coupled_model,
)
    # This is the upstream SpeedyWeather exchange sequence with one essential
    # fold-row finalization inserted between each sparse regrid and any
    # rotation/halo operation. The generic extension writes into
    # `vec(interior(field))`, bypassing ConservativeRegridding's field-aware
    # RightCenterFolded partner mirroring and leaving NaNs in duplicate slots.
    from_atmosphere = exchanger.regridder.from_atmosphere
    exchange_state = exchanger.state
    surface_layer = atmosphere.model.spectral_grid.nlayers

    u_atmosphere = SpeedyWeather.RingGrids.field_view(
        atmosphere.variables.grid.u,
        :,
        surface_layer,
    ).data
    v_atmosphere = SpeedyWeather.RingGrids.field_view(
        atmosphere.variables.grid.v,
        :,
        surface_layer,
    ).data
    temperature_atmosphere = SpeedyWeather.RingGrids.field_view(
        atmosphere.variables.grid.temperature,
        :,
        surface_layer,
    ).data
    humidity_atmosphere = SpeedyWeather.RingGrids.field_view(
        atmosphere.variables.grid.humidity,
        :,
        surface_layer,
    ).data
    pressure_atmosphere = exp.(atmosphere.variables.grid.pressure.data)
    shortwave_down =
        atmosphere.variables.parameterizations.surface_shortwave_down.data
    longwave_down =
        atmosphere.variables.parameterizations.surface_longwave_down.data
    rain = atmosphere.variables.parameterizations.rain_rate.data
    snow = haskey(atmosphere.variables.parameterizations, :snow_rate) ?
        atmosphere.variables.parameterizations.snow_rate.data : nothing

    for (destination, source) in (
        (exchange_state.u, u_atmosphere),
        (exchange_state.v, v_atmosphere),
        (exchange_state.T, temperature_atmosphere),
        (exchange_state.q, humidity_atmosphere),
        (exchange_state.p, pressure_atmosphere),
        (exchange_state.ℐꜜˢʷ, shortwave_down),
        (exchange_state.ℐꜜˡʷ, longwave_down),
        (exchange_state.Jʳⁿ, rain),
    )
        ConservativeRegridding.regrid!(
            destination,
            from_atmosphere,
            source,
        )
        _mirror_right_center_fold_partners!(destination)
    end
    if !isnothing(snow)
        ConservativeRegridding.regrid!(
            exchange_state.Jˢⁿ,
            from_atmosphere,
            snow,
        )
        _mirror_right_center_fold_partners!(exchange_state.Jˢⁿ)
    end

    arch = Oceananigans.architecture(exchange_grid)
    Oceananigans.Utils.launch!(
        arch,
        exchange_grid,
        :xy,
        _NUMERICAL_EARTH_SPEEDY_EXTENSION._rotate_winds!,
        exchange_state.u,
        exchange_state.v,
        exchange_grid,
    )

    Oceananigans.fill_halo_regions!((exchange_state.u, exchange_state.v))
    for field in (
        exchange_state.T,
        exchange_state.q,
        exchange_state.p,
        exchange_state.ℐꜜˢʷ,
        exchange_state.ℐꜜˡʷ,
        exchange_state.Jʳⁿ,
    )
        Oceananigans.fill_halo_regions!(field)
    end
    isnothing(snow) || Oceananigans.fill_halo_regions!(exchange_state.Jˢⁿ)

    _precipitation_depth_to_mass_flux!(
        exchanger.state.Jʳⁿ.data,
        exchanger.state.Jˢⁿ.data,
    )
    return nothing
end

const _ReadyESMRadiativeSpeedyEarthSystem = NumericalEarth.EarthSystemModel{
    <:AtmosphereDrivenRadiation,
    <:_ReadyESMSpeedySimulation,
}

@inline function _conditional_ocean_flux(grid_cell_mean_flux, land_fraction)
    ocean_fraction = one(land_fraction) - land_fraction
    safe_ocean_fraction = max(ocean_fraction, eps(typeof(ocean_fraction)))
    return ifelse(
        ocean_fraction > zero(ocean_fraction),
        grid_cell_mean_flux / safe_ocean_fraction,
        zero(grid_cell_mean_flux),
    )
end

"""
Copy NumericalEarth's newly regridded ocean/ice turbulent fluxes into the
prognostic input fields consumed by SpeedyWeather's prescribed-ocean surface
flux parameterizations, converting from a whole-grid-cell mean to the
conditional ocean-area mean expected by those parameterizations.

NumericalEarth currently writes these values to
`variables.parameterizations.ocean`, which are diagnostic work arrays reset at
the start of every SpeedyWeather timestep. `PrescribedOceanHeatFlux` and
`PrescribedOceanHumidityFlux` instead read `variables.prognostic.ocean`; unless
the values are synchronized here, the coupled atmosphere receives zero ocean
sensible-heat and water-vapor flux on the following timestep.

Conservative regridding has already included zero-flux inactive land cells in
the destination-cell average. SpeedyWeather subsequently multiplies its
prescribed ocean input by the destination ocean fraction, so the input must be
the conditional ocean mean to avoid applying that area fraction twice.
"""
function _synchronize_prescribed_ocean_flux_inputs!(atmosphere)
    prognostic_ocean = atmosphere.variables.prognostic.ocean
    diagnosed_ocean = atmosphere.variables.parameterizations.ocean
    land_fraction = atmosphere.model.land_sea_mask.mask.data
    prognostic_ocean.sensible_heat_flux.data .=
        _conditional_ocean_flux.(
            diagnosed_ocean.sensible_heat_flux.data,
            land_fraction,
        )
    prognostic_ocean.surface_humidity_flux.data .=
        _conditional_ocean_flux.(
            diagnosed_ocean.surface_humidity_flux.data,
            land_fraction,
        )
    return nothing
end

"""
Complete NumericalEarth's SpeedyWeather state update by preserving regridded
ocean/ice turbulent fluxes in the prescribed-flux inputs and returning dynamic
ClimaSeaIce concentration for the following atmospheric step. This keeps the
atmospheric tendencies, RRTMGP albedo, and NumericalEarth interfaces on the
same surface state.
"""
function NumericalEarth.EarthSystemModels.update_net_fluxes!(
    coupled_model::_ReadyESMRadiativeSpeedyEarthSystem,
    atmosphere::_ReadyESMSpeedySimulation,
)
    # First retain NumericalEarth's turbulent-flux and mixed ice/ocean surface-
    # temperature exchange.
    invoke(
        NumericalEarth.EarthSystemModels.update_net_fluxes!,
        Tuple{Any, _ReadyESMSpeedySimulation},
        coupled_model,
        atmosphere,
    )

    # Overwrite the generic extension's eager ocean/ice mixture with a dry-cell
    # masked, component-aware conservative map before copying the fluxes into
    # SpeedyWeather's prescribed inputs.
    _regrid_surface_fluxes!(coupled_model, atmosphere)

    # NumericalEarth's SpeedyWeather extension regrids into parameterization
    # work arrays, but the prescribed-ocean schemes read prognostic inputs on
    # the next atmospheric step. Preserve the freshly regridded mixed
    # ocean/sea-ice fluxes across that reset boundary.
    _synchronize_prescribed_ocean_flux_inputs!(atmosphere)

    _regrid_surface_intensives!(coupled_model, atmosphere)
    return nothing
end
