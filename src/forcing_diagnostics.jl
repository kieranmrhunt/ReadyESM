"""
Evaluate instantaneous clear-sky CO₂ forcing on a fixed RRTMGP standard
atmosphere. Positive `forcing_wm2` means less outgoing longwave radiation than
the lowest-CO₂ member.
"""
function run_co2_forcing_sweep(
    co2_ppm::AbstractVector{<:Real};
    output_dir::AbstractString = joinpath(PROJECT_ROOT, "artifacts", "co2_forcing_sweep"),
    nlay::Int = 60,
)
    isempty(co2_ppm) && throw(ArgumentError("co2_ppm cannot be empty"))
    all(>(0), co2_ppm) || throw(ArgumentError("all CO2 concentrations must be positive"))
    concentrations = sort(Float64.(co2_ppm))

    profile = RRTMGP.standard_atmosphere(Float64; kind = :midlatitude_summer, nlay)
    output = RRTMGP.solve(
        profile;
        method = RRTMGP.ClearSkyRadiation(false),
        cos_zenith = 0.5,
        surface_albedo = 0.2,
        surface_emissivity = 0.98,
    )
    solver = output.solver
    olr = similar(concentrations)
    surface_lw_down = similar(concentrations)
    for (i, ppm) in pairs(concentrations)
        RRTMGP.set_volume_mixing_ratio!(solver, "co2", ppm * 1.0e-6)
        RRTMGP.update_fluxes!(solver)
        olr[i] = RRTMGP.lw_flux_up(solver)[end, 1]
        surface_lw_down[i] = RRTMGP.lw_flux_dn(solver)[1, 1]
    end
    forcing = first(olr) .- olr

    mkpath(output_dir)
    netcdf_path = joinpath(output_dir, "co2_forcing_sweep.nc")
    NCDataset(netcdf_path, "c") do dataset
        defDim(dataset, "scenario", length(concentrations))
        for (name, values, units) in (
            ("co2", concentrations, "ppm"),
            ("outgoing_longwave", olr, "W/m^2"),
            ("surface_longwave_down", surface_lw_down, "W/m^2"),
            ("instantaneous_forcing", forcing, "W/m^2"),
        )
            variable = defVar(dataset, name, Float64, ("scenario",))
            variable.attrib["units"] = units
            variable[:] = values
        end
        dataset.attrib["reference_co2_ppm"] = first(concentrations)
        dataset.attrib["profile"] = "RRTMGP analytic midlatitude_summer"
        dataset.attrib["sky"] = "clear"
    end

    figure = Figure(size = (1000, 450))
    axis_olr = Axis(
        _plot_layout(figure, 1, 1);
        title = "Outgoing longwave radiation",
        xlabel = "CO₂ (ppm)",
        ylabel = "OLR (W m⁻²)",
    )
    lines!(axis_olr, concentrations, olr; linewidth = 2)
    scatter!(axis_olr, concentrations, olr; markersize = 10)
    axis_forcing = Axis(
        _plot_layout(figure, 1, 2);
        title = "Instantaneous forcing relative to $(round(Int, first(concentrations))) ppm",
        xlabel = "CO₂ (ppm)",
        ylabel = "reduction in OLR (W m⁻²)",
    )
    lines!(axis_forcing, concentrations, forcing; linewidth = 2)
    scatter!(axis_forcing, concentrations, forcing; markersize = 10)
    figure_path = joinpath(output_dir, "co2_forcing_sweep.png")
    save(figure_path, figure; px_per_unit = 1.5)

    return (;
        co2_ppm = concentrations,
        outgoing_longwave = olr,
        surface_longwave_down = surface_lw_down,
        forcing_wm2 = forcing,
        netcdf_path,
        figure_path,
    )
end
