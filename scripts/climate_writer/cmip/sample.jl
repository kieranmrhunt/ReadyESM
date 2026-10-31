# Included by the streamed runner after WriterIO, WriterSampling and the existing samplers (Claude, 2026-10-29).
# Reads only: no model state is changed. See DESIGN.md for each definition.
include("cmip_fields.jl")
CMIPFields.register_metadata!(WriterIO.WriterGeometry.SPECS)

function cmip_atmosphere_samples(R, earth)
    atmosphere = earth.atmosphere
    vars, model = atmosphere.variables, atmosphere.model
    p = vars.parameterizations
    host(field) = Array(Base.has_offset_axes(field.data) ? parent(field.data) : field.data)
    npoints = length(host(p.surface_air_temperature))
    radiation = model.longwave_radiation
    launch = earth.clock.iteration == 0
    # Before the first completed step the geopotential may hold constructor values, which need not be finite.
    # Interpolate a finite placeholder there and mark zg unavailable below (Codex review, 29 Oct).
    geopotential = host(vars.grid.geopotential)
    launch && (geopotential = zeros(size(geopotential)))
    profiles = CMIPFields.profile_samples(host(vars.grid.temperature), host(vars.grid.u), host(vars.grid.v),
        host(vars.grid.humidity), geopotential, exp.(Float64.(vec(host(vars.grid.pressure)))),
        Array(model.geometry.σ_levels_full), model.planet.gravity)
    if launch || radiation.call_counter == 0
        # Before the first completed step the parameterization arrays and the geopotential hold constructor values
        # and the cloud arrays have no solve. Initial files are outside the monthly means, so mark them unavailable.
        earth.clock.iteration == 0 || error("radiation has no solve after a completed step")
        result = Dict{String,Any}(name => fill(NaN, npoints) for (name, _, _) in CMIPFields.SURFACE_FIELDS)
        for (name, values) in profiles
            startswith(name, "cmip_zg_") && (values .= NaN)
            result[name] = values
        end
        sea = vec(host(vars.prognostic.ocean.sea_surface_temperature))
        land_t = host(vars.prognostic.land.soil_temperature)
        land_t = ndims(land_t) == 1 ? land_t : land_t[:, 1]
        land_fraction = vec(R._host_array(model.land_sea_mask.mask))
        result["cmip_ts"] = [CMIPFields.radiative_surface_temperature(land_fraction[i], sea[i], land_t[i])
                             for i in 1:npoints]
        return result
    end
    momentum = R._plain_surface_momentum_flux(model.surface_momentum_flux)
    momentum isa R.SpeedyWeather.SurfaceMomentumFlux || error("tauu needs SpeedyWeather's SurfaceMomentumFlux")
    land_fraction = vec(R._host_array(model.land_sea_mask.mask))
    drag = CMIPFields.selected_drag(momentum.use_boundary_layer_drag, vec(host(p.boundary_layer_drag)), land_fraction,
        momentum.drag_land, momentum.drag_ocean)
    u_prev, v_prev = host(vars.grid.u_prev), host(vars.grid.v_prev)
    tauu, tauv = CMIPFields.surface_stress(vec(host(p.surface_air_density)), drag, vec(host(p.surface_wind_speed)),
        momentum.wind_slowdown, u_prev[:, end], v_prev[:, end])
    land_t = host(vars.prognostic.land.soil_temperature)
    land_t = ndims(land_t) == 1 ? land_t : land_t[:, 1]
    solver = radiation.solver
    RR = R.RRTMGP
    result = Dict{String,Any}(CMIPFields.surface_samples(vec(host(p.surface_air_temperature)), land_fraction,
        vec(host(vars.prognostic.ocean.sea_surface_temperature)), land_t, Array(RR.cloud_fraction(solver)),
        Array(RR.cloud_liquid_water_path(solver)), Array(RR.cloud_ice_water_path(solver)), tauu, tauv))
    merge!(result, profiles)
    return result
end
