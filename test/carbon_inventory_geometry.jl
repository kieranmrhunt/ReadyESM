using KernelAbstractions

@kernel function _carbon_test_native_volumes!(volumes, grid)
    i,j,k = @index(Global, NTuple)
    c = ReadyESM.Oceananigans.Center()
    inactive = ReadyESM.Oceananigans.Grids.inactive_node(i,j,k,grid,c,c,c)
    @inbounds volumes[i,j,k] = inactive ? 0.0 : Float64(
        ReadyESM.Oceananigans.volume(i,j,k,grid,c,c,c))
end

function _carbon_test_volumes(grid)
    O = ReadyESM.Oceananigans
    volumes = O.on_architecture(O.Architectures.architecture(grid), zeros(Float64,size(grid)))
    backend = KernelAbstractions.get_backend(volumes)
    _carbon_test_native_volumes!(backend,64)(volumes,grid;ndrange=size(volumes))
    KernelAbstractions.synchronize(backend)
    return Array(volumes)
end

function _carbon_test_geometry_grid(arch, FT, kind, immersed)
    O = ReadyESM.Oceananigans
    z = O.Grids.MutableVerticalDiscretization((-200,0))
    grid = if kind == :rectilinear
        O.RectilinearGrid(arch,FT;size=(8,8,5),x=(0,360),y=(-75,75),z,halo=(3,3,3))
    elseif kind == :latitude_longitude
        O.LatitudeLongitudeGrid(arch,FT;size=(8,8,5),longitude=(0,360),latitude=(-75,75),z,halo=(3,3,3))
    else
        O.TripolarGrid(arch,FT;size=(20,32,5),z,halo=(5,5,4))
    end
    bottom(x,y) = x > 180 ? 20 : -120
    return immersed ? O.ImmersedBoundaryGrid(grid,O.GridFittedBottom(bottom)) : grid
end

@testset "Carbon inventory preserves evolving ocean geometry" begin
    O = ReadyESM.Oceananigans
    architectures = Any[O.CPU()]
    if get(ENV,"READYESM_GEOMETRY_DEVICE","cpu") == "gpu"
        ReadyESM.CUDA.functional(true) || error("GPU geometry tests require CUDA")
        ReadyESM.CUDA.allowscalar(false)
        push!(architectures,O.GPU())
    end
    for arch in architectures, FT in (Float32,Float64),
            kind in (:rectilinear,:latitude_longitude,:tripolar), immersed in (false,true)
        @testset "$(typeof(arch)) $FT $kind immersed=$immersed" begin
            grid = _carbon_test_geometry_grid(arch,FT,kind,immersed)
            underlying = immersed ? grid.underlying_grid : grid
            dic = O.CenterField(grid)
            model = (ocean=(model=(;grid,tracers=(DIC=dic,)),),
                interfaces=(ocean_properties=(reference_density=1024.0,),))
            initial_volumes = _carbon_test_volumes(grid)
            wet = initial_volumes .> 0
            @test any(wet)
            @test immersed ? any(.!wet) : all(wet)
            for stretch in (FT(0.93),FT(1),FT(1.07))
                for name in (:σᶜᶜⁿ,:σᶠᶜⁿ,:σᶜᶠⁿ,:σᶠᶠⁿ,:σᶜᶜ⁻)
                    fill!(getproperty(underlying.z,name),stretch)
                end
                fill!(underlying.z.ηⁿ,(stretch-one(FT))*FT(200))
                volumes = _carbon_test_volumes(grid)
                @test isapprox(sum(volumes),sum(initial_volumes)*Float64(stretch);rtol=16eps(FT))
                host_grid = ReadyESM.Adapt.adapt(Array,grid)
                host_underlying = immersed ? host_grid.underlying_grid : host_grid
                for name in (:σᶜᶜⁿ,:σᶠᶜⁿ,:σᶜᶠⁿ,:σᶠᶠⁿ,:σᶜᶜ⁻,:ηⁿ)
                    host_values = getproperty(host_underlying.z,name)
                    device_values = getproperty(underlying.z,name)
                    @test axes(host_values) == axes(device_values)
                    @test Array(parent(host_values)) == Array(parent(device_values))
                end
                for mode in (:prescribed,:tracer_only)
                    offset = ReadyESM._carbon_storage_offset(mode,Val(:DIC),Float64)
                    # Inactive cells deliberately contain NaN. The inventory
                    # must mask them before checking concentrations or summing.
                    values = fill(FT(NaN),size(grid))
                    values[wet] .= FT(2125/Float64(stretch)-offset)
                    copyto!(O.interior(dic),O.on_architecture(arch,values))
                    expected = sum((Float64.(values[wet]).+offset).*volumes[wet]) *
                        1024 * 1e-6 * ReadyESM._CARBON_MOLAR_MASS_G_MOL * ReadyESM._CARBON_PG_PER_G
                    actual = ReadyESM._carbon_ocean_inventory_pgc(model,mode)
                    @test isapprox(actual,expected;rtol=64eps(Float64),atol=0)
                    @test isequal(Array(O.interior(dic)),values)
                    @test _carbon_test_volumes(grid) == volumes
                end
            end
        end
    end
end
