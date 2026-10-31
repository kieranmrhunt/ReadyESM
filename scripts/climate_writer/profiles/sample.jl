# Include after the existing writer and surface/radiation sampler.
include("profiles.jl")
ClimateProfiles.register_metadata!(WriterIO.WriterGeometry.SPECS)

function climate_profile_samples(R, atmosphere)
    vars, model = atmosphere.variables, atmosphere.model
    host(field) = Array(Base.has_offset_axes(field.data) ? parent(field.data) : field.data)
    return ClimateProfiles.samples(host(vars.grid.temperature),host(vars.grid.humidity),host(vars.grid.u),
        exp.(Float64.(vec(host(vars.grid.pressure)))),Array(model.geometry.σ_levels_full),
        Array(model.geometry.σ_levels_thick),model.planet.gravity)
end
