"""
Speed-limit drag with an additional surface-only safeguard for externally
coupled atmosphere configurations.

`all_levels` retains SpeedyWeather's ordinary unresolved-wind protection.  The
surface branch uses the same quadratic excess-speed law, but acts only on the
lowest atmospheric layer.  It therefore leaves ordinary surface winds and the
resolved upper-level jet untouched below their respective thresholds.
"""
struct CoupledSurfaceSpeedLimitDrag{D, NF} <: SpeedyWeather.AbstractDrag
    all_levels::D
    surface_speed_limit::NF
    surface_drag::NF
end

Adapt.@adapt_structure CoupledSurfaceSpeedLimitDrag

function CoupledSurfaceSpeedLimitDrag(
    spectral_grid;
    speed_limit,
    drag,
    surface_speed_limit,
    surface_drag,
)
    all_levels = SpeedyWeather.SpeedLimitDrag(
        spectral_grid;
        speed_limit,
        drag,
    )
    NF = spectral_grid.NF
    return CoupledSurfaceSpeedLimitDrag(
        all_levels,
        convert(NF, surface_speed_limit),
        convert(NF, surface_drag),
    )
end

function SpeedyWeather.initialize!(
    scheme::CoupledSurfaceSpeedLimitDrag,
    model::SpeedyWeather.AbstractModel,
)
    SpeedyWeather.initialize!(scheme.all_levels, model)
    return nothing
end

function SpeedyWeather.drag!(
    vars,
    scheme::CoupledSurfaceSpeedLimitDrag,
    lf::Integer,
    model,
)
    SpeedyWeather.drag!(vars, scheme.all_levels, lf, model)

    surface = size(vars.grid.u, 2)
    u = SpeedyWeather.RingGrids.field_view(vars.grid.u, :, surface)
    v = SpeedyWeather.RingGrids.field_view(vars.grid.v, :, surface)
    Fu = SpeedyWeather.RingGrids.field_view(
        vars.tendencies.grid.u,
        :,
        surface,
    )
    Fv = SpeedyWeather.RingGrids.field_view(
        vars.tendencies.grid.v,
        :,
        surface,
    )

    c = scheme.surface_drag * vars.prognostic.scale[]
    speed_limit = scheme.surface_speed_limit
    @. Fu -= c * max(0, sqrt(u^2 + v^2) - speed_limit)^2 * sign(u)
    @. Fv -= c * max(0, sqrt(u^2 + v^2) - speed_limit)^2 * sign(v)
    return nothing
end
