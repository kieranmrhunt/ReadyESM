module LongRunFields

export radiation_fields, ice_fields, register_metadata!, SPECS

const SPECS = Dict(
    "ice_thickness" => (:surface, "m", "conditional sea-ice thickness"),
    "ice_volume_per_area" => (:surface, "m", "sea-ice volume per grid-cell area"),
    "toa_incoming_shortwave" => (:atmosphere, "W m-2", "incoming shortwave at TOA, latest cached solve"),
    "toa_outgoing_shortwave" => (:atmosphere, "W m-2", "outgoing shortwave at TOA, latest cached solve"),
    "toa_outgoing_longwave" => (:atmosphere, "W m-2", "outgoing longwave at TOA, latest cached solve"),
    "toa_net_downward" => (:atmosphere, "W m-2", "net downward radiation at TOA, latest cached solve"),
    "toa_clear_outgoing_shortwave" => (:atmosphere, "W m-2", "clear-sky outgoing shortwave at TOA, latest cached solve"),
    "toa_clear_outgoing_longwave" => (:atmosphere, "W m-2", "clear-sky outgoing longwave at TOA, latest cached solve"),
    "toa_clear_net_downward" => (:atmosphere, "W m-2", "clear-sky net downward radiation at TOA, latest cached solve"),
    "surface_shortwave_down" => (:atmosphere, "W m-2", "downward shortwave at surface, latest cached solve"),
    "surface_shortwave_up" => (:atmosphere, "W m-2", "upward shortwave at surface, latest cached solve"),
    "surface_longwave_down" => (:atmosphere, "W m-2", "downward longwave at surface, latest cached solve"),
    "surface_longwave_up" => (:atmosphere, "W m-2", "upward longwave at surface, latest cached solve"),
    "surface_net_radiation_downward" => (:atmosphere, "W m-2", "net downward surface radiation, latest cached solve"),
    "snow_rate" => (:atmosphere, "m s-1", "solid precipitation rate, liquid water equivalent"),
    "total_precipitation_rate" => (:atmosphere, "m s-1", "rain plus snow, liquid water equivalent"),
    "surface_sensible_heat_upward" => (:atmosphere, "W m-2", "surface sensible heat supplied to atmosphere in completed step"),
    "surface_water_vapor_upward" => (:atmosphere, "kg m-2 s-1", "surface water vapor supplied to atmosphere in completed step"),
    "global_toa_net_downward" => (:scalar, "W m-2", "global net downward TOA radiation, latest cached solve"),
    "global_total_precipitation_rate" => (:scalar, "m s-1", "global rain plus snow, liquid water equivalent"),
    "global_ice_area" => (:scalar, "m2", "ice concentration integrated over unique active ocean area"),
    "global_ice_volume" => (:scalar, "m3", "ice volume integrated over unique active ocean area"),
    "northern_ice_area" => (:scalar, "m2", "ice area north of or on the equator"),
    "southern_ice_area" => (:scalar, "m2", "ice area south of the equator"),
    "northern_ice_volume" => (:scalar, "m3", "ice volume north of or on the equator"),
    "southern_ice_volume" => (:scalar, "m3", "ice volume south of the equator"),
)

function register_metadata!(registry)
    isempty(intersect(keys(registry), keys(SPECS))) || error("extended writer names already registered")
    merge!(registry, SPECS)
    return nothing
end

"""Read RRTMGP's bottom-to-top arrays without modifying any input.
All terms come from one cached solve. They are post-step diagnostic fluxes,
not a reconstruction of the heating applied over that completed step.
"""
function radiation_fields(sw_down, sw_up, lw_down, lw_up, clear_sw_up, clear_lw_up)
    fluxes = (sw_down, sw_up, lw_down, lw_up, clear_sw_up, clear_lw_up)
    all(x -> ndims(x) == 2 && size(x) == size(sw_down), fluxes) ||
        throw(DimensionMismatch("radiation flux shapes differ"))
    size(sw_down, 1) >= 2 && size(sw_down, 2) > 0 || error("empty radiation column")
    all(x -> all(isfinite, x), fluxes) || error("nonfinite cached radiation")
    top(x) = Float64.(vec(x[end, :]))
    bottom(x) = Float64.(vec(x[1, :]))
    incoming, reflected, emitted = top(sw_down), top(sw_up), top(lw_up)
    clear_reflected, clear_emitted = top(clear_sw_up), top(clear_lw_up)
    result = Dict(
        "toa_incoming_shortwave" => incoming,
        "toa_outgoing_shortwave" => reflected,
        "toa_outgoing_longwave" => emitted,
        "toa_net_downward" => incoming .- reflected .- emitted,
        "toa_clear_outgoing_shortwave" => clear_reflected,
        "toa_clear_outgoing_longwave" => clear_emitted,
        "toa_clear_net_downward" => incoming .- clear_reflected .- clear_emitted,
        "surface_shortwave_down" => bottom(sw_down),
        "surface_shortwave_up" => bottom(sw_up),
        "surface_longwave_down" => bottom(lw_down),
        "surface_longwave_up" => bottom(lw_up),
        "surface_net_radiation_downward" => bottom(sw_down) .- bottom(sw_up) .+ bottom(lw_down) .- bottom(lw_up),
    )
    return result
end

function ice_fields(thickness, concentration)
    ndims(thickness) == 2 && size(thickness) == size(concentration) ||
        throw(DimensionMismatch("ice field shapes differ"))
    # Multiply in the native precision before widening, as in the mass ledger.
    return Dict("ice_thickness" => Float64.(thickness),
        "ice_volume_per_area" => Float64.(thickness .* concentration))
end

end
