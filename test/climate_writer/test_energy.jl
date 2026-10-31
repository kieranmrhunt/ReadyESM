using Test
const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(PROJECT_ROOT, "scripts", "climate_writer", "energy", "sample.jl"))
const E = ClimateEnergy

struct GuardedArray{N} <: AbstractArray{Float64,N}
    values::Array{Float64,N}
    allowed::AbstractArray{Bool,N}
end
Base.size(x::GuardedArray) = size(x.values)
Base.IndexStyle(::Type{<:GuardedArray}) = IndexLinear()
Base.getindex(x::GuardedArray, i::Int) = x.allowed[i] ? x.values[i] : error("absent tile read")

@testset "Native latent heat components" begin
    mask = reshape([1,1,1,0,0,1],3,2)
    a = reshape([0.,1.,.25,NaN,NaN,.5],3,2)
    water = reshape([100.,NaN,80.,NaN,NaN,-20.],3,2)
    ice = reshape([NaN,60.,40.,NaN,NaN,-10.],3,2)
    old = deepcopy((mask,a,water,ice))
    guarded_water = GuardedArray(water, (mask .== 1) .& (a .< 1))
    guarded_ice = GuardedArray(ice, (mask .== 1) .& (a .> 0))
    guarded_a = GuardedArray(a, mask .== 1)
    ocean = E.ocean_samples(guarded_water, guarded_ice, guarded_a, mask)
    @test ocean["ocean_ice_latent_heat_upward"] == reshape([100.,60.,70.,0.,0.,-15.],3,2)
    @test ocean["open_water_latent_heat_contribution"] == reshape([100.,0.,60.,0.,0.,-10.],3,2)
    @test ocean["sea_ice_latent_heat_contribution"] == reshape([0.,60.,10.,0.,0.,-5.],3,2)
    @test ocean["open_water_fraction"] + ocean["sea_ice_surface_fraction"] == mask
    vapour = E.vapor_samples(guarded_water,guarded_ice,guarded_a,mask)
    @test vapour["open_water_vapor_contribution"] == ocean["open_water_latent_heat_contribution"]
    @test vapour["sea_ice_vapor_contribution"] == ocean["sea_ice_latent_heat_contribution"]
    @test isequal(old, (mask,a,water,ice))
    for invalid in (-.01,1.01,NaN,Inf)
        bad = copy(a)
        bad[1] = invalid
        @test_throws ErrorException E.ocean_samples(water,ice,bad,mask)
    end
    for invalid in (NaN,Inf,-Inf)
        bad = copy(water)
        bad[1] = invalid
        @test_throws ErrorException E.ocean_samples(bad,ice,a,mask)
        bad = copy(ice)
        bad[2] = invalid
        @test_throws ErrorException E.ocean_samples(water,bad,a,mask)
    end
    @test_throws ErrorException E.ocean_samples(water,ice,a,fill(.5,3,2))
    @test_throws DimensionMismatch E.ocean_samples(water[1:2,:],ice,a,mask)
    @test_throws ErrorException E.ocean_samples(zeros(0,0),zeros(0,0),zeros(0,0),zeros(0,0))

    q = [1e-4,NaN,-2e-5,0.]
    fraction = [.5,0.,1.,.25]
    latent = 2.5e6
    land = E.land_samples(GuardedArray(q,fraction .> 0),fraction,latent)
    @test land["land_latent_heat_contribution"] ≈ [125.,0.,-50.,0.]
    @test land["land_surface_fraction"] == fraction
    @test_throws ErrorException E.land_samples(q,fraction,latent;snow_enabled=true)
    for bad in (0.,-1.,NaN,Inf)
        @test_throws ErrorException E.land_samples(q,fraction,bad)
    end
    for bad in (-.1,1.1,NaN,Inf)
        @test_throws ErrorException E.land_samples(q,[bad,0.,1.,.25],latent)
    end
    @test_throws ErrorException E.land_samples([NaN,NaN,0.,0.],fraction,latent)
    @test_throws DimensionMismatch E.land_samples(q[1:3],fraction,latent)
    @test q[1] == 1e-4 && fraction == [.5,0.,1.,.25]
    for F in (Float32,Float64), heat in (-100.,0.,80.,400.)
        stored = F(heat/latent)
        reconstructed = E.land_samples([stored],[1.],latent)["land_latent_heat_contribution"][1]
        @test isapprox(reconstructed,heat;rtol=2eps(F),atol=0.)
    end

    area = ones(3,2)
    w = fill(.25,4)
    totals = E.global_summary(ocean,land,area,w,16.)
    @test totals.ocean_power_w == 215.
    @test totals.land_power_w ≈ 300.
    @test totals.global_latent_heat_upward ≈ 515/16
    @test totals.represented_surface_fraction == 4/16 + 1.75/4
    @test totals.represented_surface_fraction != 1
    @test_throws ErrorException E.global_summary(ocean,land,area,w,0.)
    @test_throws ErrorException E.global_summary(ocean,land,area,ones(4),16.)
    @test_throws DimensionMismatch E.global_summary(ocean,land,ones(2,2),w,16.)
    for F in (Float32,Float64)
        o = E.ocean_samples(fill(F(100),2,2),fill(F(80),2,2),fill(F(.25),2,2),ones(2,2))
        @test all(==(95.),o["ocean_ice_latent_heat_upward"])
    end
    specs = Dict{String,Tuple{Symbol,String,String}}()
    E.register_metadata!(specs)
    @test Set(keys(specs)) == Set(keys(ocean)) ∪ Set(keys(land)) ∪ Set(keys(vapour))
    @test_throws ErrorException E.register_metadata!(specs)
end

@testset "Read-only adapter and availability" begin
    mock_R = (; _host_array=copy, _interior_matrix=copy)
    land_model = (; snow=nothing, constants=(;thermodynamics=(;latent_heat_vaporization=2.5e6)))
    atmosphere = (; model=(;land=(;model=land_model),land_sea_mask=(;mask=[1.,0.])),
        variables=(;prognostic=(;land=(;surface_humidity_flux=[1e-4,NaN]))))
    interfaces = (; atmosphere_ocean_interface=(;fluxes=(;latent_heat=fill(100.,2,1),water_vapor=fill(1e-5,2,1))),
        atmosphere_sea_ice_interface=(;fluxes=(;latent_heat=fill(80.,2,1),water_vapor=fill(-2e-6,2,1))))
    earth = (; atmosphere, interfaces, clock=(;iteration=1),
        sea_ice=(;model=(;ice_concentration=reshape([0.,.5],2,1))))
    before = deepcopy(earth)
    fields = climate_energy_samples(mock_R,earth,ones(2,1))
    @test fields["land_latent_heat_contribution"] == [250.,0.]
    @test fields["ocean_ice_latent_heat_upward"] == reshape([100.,90.],2,1)
    @test fields["open_water_vapor_contribution"] == reshape([1e-5,5e-6],2,1)
    @test fields["sea_ice_vapor_contribution"] == reshape([0.,-1e-6],2,1)
    @test isequal(before,earth)
    initial = climate_energy_samples(mock_R,merge(earth,(;clock=(;iteration=0))),ones(2,1))
    @test Set(keys(initial)) == Set(keys(fields))
    @test all(x -> all(isnan,x),values(initial))
    # Constructor fluxes must never be read on the unavailable launch path.
    minimal = (; atmosphere, clock=(;iteration=0))
    @test all(x -> all(isnan,x),values(climate_energy_samples(mock_R,minimal,ones(2,1))))
end
println("CLIMATE_ENERGY_COMPONENT_TEST_PASS")
