#!/usr/bin/env julia

using CairoMakie
using NCDatasets
using Printf

# Chart contract:
# - Question: what is the coupled model's actual final saved surface state,
#   and do the principal ocean, ice, atmosphere, and land fields expose a
#   localized scientific failure that a global mean would hide?
# - Takeaway target: show nine final-state maps on explicit scales, label the
#   actual saved day and scientific-audit status, and report key final values.
# - Form: a three-by-three atlas. A separate progress figure carries the full
#   trajectory; these maps answer the final-state question without compressing
#   two years of history into the same panels.
# - Data sufficiency: one complete global snapshot on the native 360x180
#   curvilinear ocean grid and T31 reduced Gaussian atmosphere grid. If the run
#   ends early, the title reports the true last saved day rather than day 730.
# - Palette: two non-neutral roots (blue and orange) plus grey inactive cells;
#   physical extrema are written into titles so interpretation does not depend
#   on colour alone. Zero ice is white and nonzero ice starts in pale blue.
# - Surface: reproducible static PNG, followed by inspection of the exported
#   image. The NetCDF path remains the canonical quantitative source.

const SECONDS_PER_DAY = 86_400.0
const INACTIVE_COLOR = colorant"#9CA3AF"
const TEMPERATURE_COLORMAP = [
    colorant"#174F7A",
    colorant"#5B91BD",
    colorant"#D7E8F3",
    colorant"#F7F7F5",
    colorant"#F6C988",
    colorant"#D97706",
    colorant"#8A4304",
]
const BLUE_COLORMAP = [
    colorant"#FFFFFF",
    colorant"#DCEAF5",
    colorant"#9ECAE1",
    colorant"#4292C6",
    colorant"#08519C",
    colorant"#08306B",
]
const ORANGE_COLORMAP = [
    colorant"#FFF7ED",
    colorant"#FED7AA",
    colorant"#FB923C",
    colorant"#C2410C",
    colorant"#7C2D12",
]

2 <= length(ARGS) <= 3 || error(
    "usage: plot_final_day_summary.jl DIAGNOSTICS.nc OUTPUT.png [RUN_STATUS.txt]",
)

diagnostics_path = abspath(ARGS[1])
output_path = abspath(ARGS[2])
status_path = length(ARGS) == 3 ? abspath(ARGS[3]) : ""
isfile(diagnostics_path) || error("diagnostics file does not exist: $diagnostics_path")

function read_variable(dataset, name)
    haskey(dataset, name) || error("diagnostics lack required variable $name")
    variable = dataset[name]
    indices = ntuple(_ -> Colon(), ndims(variable))
    return Float64.(Array(variable[indices...]))
end

function read_raw_variable(dataset, name)
    haskey(dataset, name) || error("diagnostics lack required variable $name")
    variable = NCDatasets.variable(dataset, name)
    indices = ntuple(_ -> Colon(), ndims(variable))
    return Float64.(Array(variable[indices...]))
end

diagnostics = NCDataset(diagnostics_path, "r") do dataset
    return (
        time = read_raw_variable(dataset, "time"),
        atmosphere_time =
            read_raw_variable(dataset, "atmosphere_diagnostic_time"),
        coupled_time = read_raw_variable(dataset, "coupled_budget_time"),
        ocean_longitude = read_variable(dataset, "ocean_longitude"),
        ocean_latitude = read_variable(dataset, "ocean_latitude"),
        ocean_active = read_variable(dataset, "ocean_surface_active_mask"),
        ocean_temperature = read_variable(dataset, "ocean_surface_temperature"),
        ocean_salinity = read_variable(dataset, "ocean_surface_salinity"),
        sea_ice_concentration = read_variable(dataset, "sea_ice_concentration"),
        sea_ice_thickness = read_variable(dataset, "sea_ice_thickness"),
        atmosphere_longitude = read_variable(dataset, "atmosphere_longitude"),
        atmosphere_latitude = read_variable(dataset, "atmosphere_latitude"),
        atmosphere_surface_temperature =
            read_variable(dataset, "atmosphere_surface_temperature"),
        atmosphere_rainfall_flux =
            read_variable(dataset, "atmosphere_rainfall_flux"),
        atmosphere_snowfall_flux =
            read_variable(dataset, "atmosphere_snowfall_flux"),
        column_cloud_fraction = read_variable(dataset, "column_cloud_fraction"),
        outgoing_longwave = read_variable(dataset, "outgoing_longwave"),
        land_surface_moisture = read_variable(dataset, "land_surface_moisture"),
        toa_net_downward = read_variable(dataset, "toa_net_downward"),
        global_surface_air_temperature =
            read_variable(dataset, "global_surface_air_temperature"),
        global_rainfall_flux = read_variable(dataset, "global_rainfall_flux"),
        global_snowfall_flux = read_variable(dataset, "global_snowfall_flux"),
        ocean_temperature_maximum =
            read_variable(dataset, "ocean_temperature_maximum_series"),
        ocean_salinity_minimum =
            read_variable(dataset, "ocean_salinity_minimum_series"),
        sea_ice_area_equivalent_maximum =
            read_variable(dataset, "sea_ice_max_area_equivalent_thickness"),
        land_water_budget_residual =
            read_variable(dataset, "global_land_water_budget_residual"),
    )
end

all(!isempty, (diagnostics.time, diagnostics.atmosphere_time, diagnostics.coupled_time)) ||
    error("diagnostic time axes must be nonempty")
all(issorted, (diagnostics.time, diagnostics.atmosphere_time, diagnostics.coupled_time)) ||
    error("diagnostic time axes must be monotonic")
final_day = last(diagnostics.coupled_time)
isapprox(last(diagnostics.time), final_day; atol = 1e-6, rtol = 0) || error(
    "radiation and coupled diagnostics end at different days",
)
isapprox(last(diagnostics.atmosphere_time), final_day; atol = 1e-6, rtol = 0) ||
    error("atmosphere and coupled diagnostics end at different days")

ocean_shape = size(diagnostics.ocean_temperature)
all(size(values) == ocean_shape for values in (
    diagnostics.ocean_salinity,
    diagnostics.ocean_active,
    diagnostics.sea_ice_concentration,
    diagnostics.sea_ice_thickness,
)) || error("ocean final-state fields do not share one horizontal grid")
all(value -> value == 0 || value == 1, diagnostics.ocean_active) ||
    error("ocean surface active mask is not binary")
active_ocean = diagnostics.ocean_active .== 1
any(active_ocean) || error("ocean surface active mask contains no wet cells")

function coordinate_matrices(longitude, latitude, shape)
    if ndims(longitude) == 2
        size(longitude) == shape || error("curvilinear longitude shape differs from data")
        size(latitude) == shape || error("curvilinear latitude shape differs from data")
        return longitude, latitude
    end
    length(longitude) == shape[1] || error("longitude dimension differs from data")
    length(latitude) == shape[2] || error("latitude dimension differs from data")
    return (
        repeat(reshape(longitude, :, 1), 1, length(latitude)),
        repeat(reshape(latitude, 1, :), length(longitude), 1),
    )
end

ocean_longitude, ocean_latitude = coordinate_matrices(
    diagnostics.ocean_longitude,
    diagnostics.ocean_latitude,
    ocean_shape,
)
ocean_longitude = mod.(ocean_longitude, 360)

function finite_extrema(values, mask = trues(size(values)))
    selected = values[mask .& isfinite.(values)]
    isempty(selected) && error("plot field contains no finite selected values")
    minimum_value, maximum_value = extrema(selected)
    if minimum_value == maximum_value
        padding = max(abs(minimum_value) * 0.01, 1e-6)
        return (minimum_value - padding, maximum_value + padding)
    end
    return (minimum_value, maximum_value)
end

function parse_status(path)
    isempty(path) && return "not supplied"
    isfile(path) || return "not yet audited"
    values = Dict{String, String}()
    for line in eachline(path)
        occursin('=', line) || continue
        key, value = split(line, '='; limit = 2)
        values[strip(key)] = strip(value)
    end
    return get(values, "scientific_acceptance", "not recorded")
end

function map_axis(figure, position, title; backgroundcolor = INACTIVE_COLOR)
    return Axis(
        figure[position...];
        title,
        xlabel = "longitude (°E)",
        ylabel = "latitude (°N)",
        limits = ((0, 360), (-90, 90)),
        xticks = 0:60:360,
        yticks = -60:30:60,
        backgroundcolor,
    )
end

function ocean_map!(
    figure,
    axis_position,
    colorbar_position,
    values,
    title,
    colorbar_label;
    colormap,
    colorrange = finite_extrema(values, active_ocean),
)
    axis = map_axis(figure, axis_position, title)
    shown = active_ocean .& isfinite.(values)
    plot = scatter!(
        axis,
        vec(ocean_longitude[shown]),
        vec(ocean_latitude[shown]);
        marker = :rect,
        markersize = 3.5,
        color = vec(values[shown]),
        colormap,
        colorrange,
    )
    Colorbar(figure[colorbar_position...], plot; label = colorbar_label)
    return axis
end

function atmosphere_map!(
    figure,
    axis_position,
    colorbar_position,
    values,
    title,
    colorbar_label;
    colormap,
    colorrange = finite_extrema(values),
    selected = isfinite.(values),
    colorscale = identity,
)
    axis = map_axis(figure, axis_position, title; backgroundcolor = colorant"#E5E7EB")
    shown = selected .& isfinite.(values)
    any(shown) || error("atmosphere plot field contains no finite selected values")
    plot = scatter!(
        axis,
        mod.(diagnostics.atmosphere_longitude[shown], 360),
        diagnostics.atmosphere_latitude[shown];
        marker = :rect,
        markersize = 7,
        color = values[shown],
        colormap,
        colorrange,
        colorscale,
    )
    Colorbar(figure[colorbar_position...], plot; label = colorbar_label)
    return axis
end

ocean_temperature_range =
    finite_extrema(diagnostics.ocean_temperature, active_ocean)
atmosphere_temperature_range =
    finite_extrema(diagnostics.atmosphere_surface_temperature .- 273.15)
salinity_range = finite_extrema(diagnostics.ocean_salinity, active_ocean)
area_equivalent_ice =
    diagnostics.sea_ice_concentration .* diagnostics.sea_ice_thickness
ice_maximum = maximum(area_equivalent_ice[active_ocean])
ice_range = (0.0, max(ice_maximum, 0.1))
precipitation = SECONDS_PER_DAY .* max.(
    diagnostics.atmosphere_rainfall_flux .+
    diagnostics.atmosphere_snowfall_flux,
    0,
)
precipitation_floor = 0.01
precipitation_for_plot = max.(precipitation, precipitation_floor)
precipitation_range = (
    precipitation_floor,
    max(maximum(precipitation_for_plot), precipitation_floor * 10),
)
land = isfinite.(diagnostics.land_surface_moisture)

expected_days = parse(Float64, get(ENV, "READYESM_EXPECTED_DAYS", "730"))
expected_days > 0 || error("READYESM_EXPECTED_DAYS must be positive")
completion_label = isapprox(final_day, expected_days; atol = 1e-6, rtol = 0) ?
    "complete requested horizon" :
    @sprintf("incomplete: expected day %.0f", expected_days)
acceptance = parse_status(status_path)
global_precipitation = SECONDS_PER_DAY * (
    last(diagnostics.global_rainfall_flux) +
    last(diagnostics.global_snowfall_flux)
)

figure = Figure(size = (2250, 1800), backgroundcolor = :white)
Label(
    figure[1, 1:6],
    @sprintf("ReadyESM coupled final-state summary — saved day %.2f", final_day);
    fontsize = 30,
    font = :bold,
    halign = :left,
)
Label(
    figure[2, 1:6],
    @sprintf(
        "%s · scientific acceptance: %s · TOA net %+0.2f W m⁻² · global surface air %.2f K · global precipitation %.2f mm day⁻¹",
        completion_label,
        acceptance,
        last(diagnostics.toa_net_downward),
        last(diagnostics.global_surface_air_temperature),
        global_precipitation,
    );
    fontsize = 17,
    color = colorant"#374151",
    halign = :left,
    tellwidth = false,
)

ocean_map!(
    figure, (3, 1), (3, 2), diagnostics.ocean_temperature,
    @sprintf("A  Sea-surface temperature (max %.2f °C)", maximum(diagnostics.ocean_temperature[active_ocean])),
    "°C";
    colormap = TEMPERATURE_COLORMAP,
    colorrange = ocean_temperature_range,
)
salinity_axis = ocean_map!(
    figure, (3, 3), (3, 4), diagnostics.ocean_salinity,
    @sprintf("B  Sea-surface salinity (min %.2f psu)", minimum(diagnostics.ocean_salinity[active_ocean])),
    "psu";
    colormap = BLUE_COLORMAP,
    colorrange = salinity_range,
)
active_indices = findall(active_ocean)
salinity_minimum_index = active_indices[argmin(
    diagnostics.ocean_salinity[active_ocean],
)]
scatter!(
    salinity_axis,
    [ocean_longitude[salinity_minimum_index]],
    [ocean_latitude[salinity_minimum_index]];
    marker = :xcross,
    markersize = 16,
    strokewidth = 3,
    color = colorant"#111827",
)
ocean_map!(
    figure, (3, 5), (3, 6), diagnostics.sea_ice_concentration,
    "C  Sea-ice concentration", "fraction";
    colormap = BLUE_COLORMAP,
    colorrange = (0, 1),
)

ice_axis = ocean_map!(
    figure, (4, 1), (4, 2), area_equivalent_ice,
    @sprintf("D  Area-equivalent sea-ice thickness (max %.2f m)", ice_maximum),
    "A × h (m)";
    colormap = BLUE_COLORMAP,
    colorrange = ice_range,
)
ice_maximum_index = active_indices[argmax(area_equivalent_ice[active_ocean])]
scatter!(
    ice_axis,
    [ocean_longitude[ice_maximum_index]],
    [ocean_latitude[ice_maximum_index]];
    marker = :xcross,
    markersize = 16,
    strokewidth = 3,
    color = colorant"#111827",
)
atmosphere_map!(
    figure, (4, 3), (4, 4),
    diagnostics.atmosphere_surface_temperature .- 273.15,
    "E  Surface-air temperature", "°C";
    colormap = TEMPERATURE_COLORMAP,
    colorrange = atmosphere_temperature_range,
)
atmosphere_map!(
    figure, (4, 5), (4, 6), precipitation_for_plot,
    @sprintf("F  Instantaneous precipitation (max %.1f mm day⁻¹)", maximum(precipitation)),
    "mm day⁻¹ (log colour)";
    colormap = BLUE_COLORMAP,
    colorrange = precipitation_range,
    colorscale = log10,
)

atmosphere_map!(
    figure, (5, 1), (5, 2), diagnostics.column_cloud_fraction,
    "G  Column cloud fraction", "fraction";
    colormap = BLUE_COLORMAP,
    colorrange = (0, 1),
)
atmosphere_map!(
    figure, (5, 3), (5, 4), diagnostics.land_surface_moisture,
    "H  Land surface soil saturation", "fraction";
    colormap = BLUE_COLORMAP,
    colorrange = (0, 1),
    selected = land,
)
atmosphere_map!(
    figure, (5, 5), (5, 6), diagnostics.outgoing_longwave,
    "I  Outgoing longwave radiation", "W m⁻²";
    colormap = ORANGE_COLORMAP,
)

Label(
    figure[6, 1:6],
    @sprintf(
        "Final full-depth ocean Tmax %.2f °C · Smin %.2f psu · max physical A×h %.2f m · land water closure residual %+0.4g kg m⁻².\nPrecipitation is the saved instantaneous flux expressed per day; grey denotes inactive/unplotted cells; crosses mark the surface salinity minimum and A×h maximum. Source: %s",
        last(diagnostics.ocean_temperature_maximum),
        last(diagnostics.ocean_salinity_minimum),
        last(diagnostics.sea_ice_area_equivalent_maximum),
        last(diagnostics.land_water_budget_residual),
        relpath(diagnostics_path, pwd()),
    );
    fontsize = 14,
    color = colorant"#4B5563",
    halign = :left,
    tellwidth = false,
)

for column in (1, 3, 5)
    colsize!(figure.layout, column, Relative(0.31))
end
for column in (2, 4, 6)
    colsize!(figure.layout, column, Relative(0.023))
end
rowgap!(figure.layout, 12)

mkpath(dirname(output_path))
save(output_path, figure; px_per_unit = 1.25)
println(
    @sprintf(
        "FINAL_DAY_SUMMARY_PASS day=%.6f expected_days=%.6f acceptance=%s output=%s",
        final_day,
        expected_days,
        acceptance,
        output_path,
    ),
)
