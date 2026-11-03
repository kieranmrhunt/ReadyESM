# Opt-in sea-ice lead-closing thickness (design artifacts/development_wave_20261103/sea_ice_lead_closing_claude/DESIGN.md).
#
# ClimaSeaIce's ProportionalEvolution grows concentration at (1 - A)/h_c times a positive cell-net volume tendency,
# with h_c the ice consolidation thickness (0.05 m here), and removes it at A/(2h) times a negative one. Hibler (1979)
# uses a thick/thin demarcation thickness h_0 of 0.5 m in the Arctic, and the literature uses 0.1-1.0 m. With a small
# h_c, alternating sub-daily freeze and melt steps can raise A under net melt (a ratchet). The option keeps the rule and
# only gives its freezing term a separate thickness h_0 (vendored ClimaSeaIce LeadClosingEvolution). It does not change
# which growth sources count, so it is a sensitivity of h_0, not a full Hibler partition.
const SEA_ICE_LEAD_CLOSING_RANGE_M = (0.05, 2.0)
const SEA_ICE_CONSOLIDATION_THICKNESS_M = 0.05

"LeadClosingEvolution for a positive setting, or `nothing` for the release rule and constructor."
function _sea_ice_concentration_evolution(config, ::Type{FT}) where {FT}
    h0 = config.sea_ice_lead_closing_thickness_m
    iszero(h0) && return nothing
    return ClimaSeaIce.SeaIceThermodynamics.LeadClosingEvolution(FT(h0))
end

"Provenance of the setting, added only when it is used, so default output is unchanged."
function _sea_ice_lead_closing_provenance(config)
    h0 = config.sea_ice_lead_closing_thickness_m
    iszero(h0) && return (;)
    return (;
        sea_ice_concentration_evolution = "proportional_with_lead_closing_thickness_v1",
        sea_ice_lead_closing_thickness_m = h0,
        sea_ice_consolidation_thickness_m = SEA_ICE_CONSOLIDATION_THICKNESS_M,
    )
end
