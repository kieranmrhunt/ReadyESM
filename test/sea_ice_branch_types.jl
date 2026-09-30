module SeaIceBranchTypesTest

using ReadyESM, Test
const O = ReadyESM.Oceananigans
const C = ReadyESM.ClimaSeaIce

function run_case(FT, DT, ice)
    grid = O.RectilinearGrid(O.CPU(), FT;
        size=(8, 8, 1), x=(0, 80_000), y=(0, 80_000), z=(-10, 0),
        halo=(4, 4, 4), topology=(O.Periodic, O.Periodic, O.Bounded))
    free_drift = (u=O.Fields.ConstantField(FT(0.07)),
                  v=O.Fields.ConstantField(FT(-0.04)))
    dynamics = C.SeaIceMomentumEquation(grid;
        rheology=C.ElastoViscoPlasticRheology(FT),
        solver=C.SplitExplicitSolver(grid; substeps=10), free_drift)
    model = C.SeaIceModel(grid; dynamics, timestepper=:ForwardEuler,
                         ice_thermodynamics=nothing, advection=nothing)
    thickness, concentration = ice === :active ? (FT(1), FT(0.9)) :
                               ice === :marginal ? (FT(0.01), FT(1e-4)) :
                                                  (zero(FT), zero(FT))
    u = reshape(FT[0.02sin(2pi*(i-1)/8)*cos(2pi*(j-1)/8) for i in 1:8, j in 1:8],8,8,1)
    v = reshape(FT[-0.01cos(2pi*(i-1)/8)*sin(2pi*(j-1)/8) for i in 1:8, j in 1:8],8,8,1)
    O.set!(model; h=thickness, ℵ=concentration, u, v)
    O.TimeSteppers.update_state!(model)
    result = Dict{String,Any}()
    for step in 1:3
        C.SeaIceDynamics.time_step_momentum!(model, model.dynamics, DT(60))
        fields = merge(model.velocities, model.dynamics.auxiliaries.fields,
                       (h=model.ice_thickness, concentration=model.ice_concentration))
        for name in (:u, :v, :σ₁₁, :σ₂₂, :σ₁₂, :h, :concentration)
            values = copy(Array(O.interior(getproperty(fields,name))))
            @test all(isfinite, values)
            result["$(step)_$(name)"] = values
        end
        if ice === :marginal
            @test all(==(FT(0.07)), result["$(step)_u"])
            @test all(==(FT(-0.04)), result["$(step)_v"])
        elseif ice === :none
            @test all(iszero, result["$(step)_u"])
            @test all(iszero, result["$(step)_v"])
        else
            @test maximum(abs, result["$(step)_σ₁₁"]) > 0
            @test maximum(abs, result["$(step)_u"]) > 0
        end
        @test all(==(thickness), result["$(step)_h"])
        @test all(==(concentration), result["$(step)_concentration"])
    end
    return result
end

cases = [(FT,DT,ice) for (FT,DT) in
         ((Float32,Float32),(Float32,Float64),(Float64,Float64))
         for ice in (:active,:marginal,:none)]

@testset "Sea-ice momentum branch types" begin
    @test startswith(realpath(pathof(C)), realpath(joinpath(ReadyESM.PROJECT_ROOT, "vendor", "ClimaSeaIce")))
    for case in cases
        run_case(case...)
    end
end

end
