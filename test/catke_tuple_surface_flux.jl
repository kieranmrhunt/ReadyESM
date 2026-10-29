module CATKETupleSurfaceFluxTests

using Test
using ReadyESM
using Oceananigans
import KernelAbstractions
using KernelAbstractions: @kernel, @index
using Oceananigans.BoundaryConditions: getbc
using Oceananigans.Models: fields
import Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities as TKE

@kernel function boundary_flux!(out, condition, grid, clock, model_fields, closures, buoyancy)
    i, j = @index(Global, NTuple)
    @inbounds out[i, j] = getbc(condition, i, j, grid, clock, model_fields, closures, buoyancy)
end

function evaluate_boundary(model, closures)
    arch = Oceananigans.Architectures.architecture(model.grid)
    output = Oceananigans.Architectures.on_architecture(arch, zeros(4, 4))
    backend = KernelAbstractions.get_backend(output)
    boundary_flux!(backend, (4, 4))(
        output, model.tracers.e.boundary_conditions.top, model.grid,
        model.clock, fields(model), closures, model.buoyancy;
        ndrange=(4, 4),
    )
    KernelAbstractions.synchronize(backend)
    return Array(output)
end

function run_tests(arch)
    @testset "Native CATKE surface flux with production closure tuples" begin
        grid = LatitudeLongitudeGrid(arch; size=(4, 4, 4), halo=(3, 3, 3),
                                     longitude=(0, 360), latitude=(-60, 60), z=(-40, 0))
        base = ReadyESM.ClimaOcean.OceanConfigurations.default_one_degree_closure()
        @test length(base) == 4
        kv, kh = CenterField(grid), CenterField(grid)
        set!(kv, 0)
        set!(kh, 0)
        extra = ReadyESM.LocalizedEstuaryDiffusivity(kv, kh)
        production = (base, (base..., extra))
        for closures in production, (stress, buoyancy_flux) in
            ((1e-4, 0.0), (0.0, 1e-7), (1e-4, 1e-7), (0.0, 0.0))
            model = HydrostaticFreeSurfaceModel(grid;
                closure=closures, buoyancy=BuoyancyTracer(), tracers=:b,
                boundary_conditions=(
                    u=FieldBoundaryConditions(top=FluxBoundaryCondition(stress)),
                    b=FieldBoundaryConditions(top=FluxBoundaryCondition(buoyancy_flux)),
                ),
            )
            # Use the actual closure tuple after constructor regularization.
            native = model.closure
            @test length(native) == length(closures)
            catke = native[1]
            coefficients = catke.turbulent_kinetic_energy_equation
            expected = -coefficients.Cᵂu★ * abs(stress)^1.5 -
                       coefficients.CᵂwΔ * max(buoyancy_flux, 0) * 10
            @test all(isapprox.(evaluate_boundary(model, native), expected; rtol=1e-13, atol=0))
            @test evaluate_boundary(model, native) == evaluate_boundary(model, catke)
            if stress > 0 || buoyancy_flux > 0
                @test expected < 0
            else
                @test expected == 0
            end
            # Existing upstream short-tuple methods remain selected.
            parameters = model.tracers.e.boundary_conditions.top.condition.parameters
            for n in 1:3
                short = (catke, ntuple(_ -> nothing, n - 1)...)
                method = which(TKE.top_tke_flux, Tuple{
                    Int, Int, typeof(grid), typeof(model.clock), typeof(fields(model)),
                    typeof(parameters), typeof(short), typeof(model.buoyancy),
                })
                @test method.module === TKE
                @test evaluate_boundary(model, short) == evaluate_boundary(model, catke)
            end
            # Position must not matter, and non-TKE closures must add zero.
            for position in eachindex(native)
                moved = (Base.tail(native)[1:position-1]..., catke,
                         Base.tail(native)[position:end]...)
                @test evaluate_boundary(model, moved) == evaluate_boundary(model, catke)
            end
            @test all(iszero, evaluate_boundary(model, (Base.tail(native)..., nothing)))
            @test all(iszero, evaluate_boundary(model, ()))
        end
    end
    println("CATKE_TUPLE_SURFACE_FLUX_PASS")
end

device = get(ENV, "READYESM_TEST_DEVICE", "cpu")
device in ("cpu", "gpu") || error("Unknown test device: $device")
if device == "gpu"
    ReadyESM.CUDA.functional() || error("CUDA unavailable")
    ReadyESM.CUDA.allowscalar(false)
end
run_tests(device == "gpu" ? GPU() : CPU())

end
