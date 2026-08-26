#!/usr/bin/env julia

using Dates
using NCDatasets
using SHA

length(ARGS) == 3 || error(
    "usage: prepare_ecco_monthly_bundle.jl RAW_DIRECTORY OUTPUT_DIRECTORY YYYY-MM",
)

raw_directory = abspath(ARGS[1])
output_directory = abspath(ARGS[2])
month_match = match(r"^(\d{4})-(\d{2})$", ARGS[3])
isnothing(month_match) && error("month must have form YYYY-MM")
source_month = Date(
    parse(Int, month_match.captures[1]),
    parse(Int, month_match.captures[2]),
    1,
)
next_month = source_month + Month(1)
month_label = Dates.format(source_month, dateformat"yyyy-mm")
filename_label = replace(month_label, "-" => "_")

isdir(raw_directory) || error("raw ECCO directory does not exist: $raw_directory")

collections = (
    temperature_salinity = (
        id = "ECCO_L4_TEMP_SALINITY_05DEG_MONTHLY_V4R4",
        filename = "OCEAN_TEMPERATURE_SALINITY_mon_mean_$(month_label)_ECCO_V4r4_latlon_0p50deg.nc",
        variables = ("THETA", "SALT"),
    ),
    velocity = (
        id = "ECCO_L4_OCEAN_VEL_05DEG_MONTHLY_V4R4",
        filename = "OCEAN_VELOCITY_mon_mean_$(month_label)_ECCO_V4r4_latlon_0p50deg.nc",
        variables = ("EVEL", "NVEL"),
    ),
    free_surface = (
        id = "ECCO_L4_SSH_05DEG_MONTHLY_V4R4",
        filename = "SEA_SURFACE_HEIGHT_mon_mean_$(month_label)_ECCO_V4r4_latlon_0p50deg.nc",
        variables = ("SSH",),
    ),
    sea_ice = (
        id = "ECCO_L4_SEA_ICE_CONC_THICKNESS_05DEG_MONTHLY_V4R4",
        filename = "SEA_ICE_CONC_THICKNESS_mon_mean_$(month_label)_ECCO_V4r4_latlon_0p50deg.nc",
        variables = ("SIheff", "SIarea"),
    ),
)

function sha512_hex(path)
    return open(path) do io
        bytes2hex(sha512(io))
    end
end

function published_sha512(path)
    checksum_path = path * ".sha512"
    isfile(checksum_path) || error("published SHA-512 file is missing: $checksum_path")
    fields = split(strip(read(checksum_path, String)))
    length(fields) >= 2 || error("malformed SHA-512 file: $checksum_path")
    lowercase(fields[1]) == fields[1] || error(
        "published SHA-512 must use lowercase hexadecimal: $checksum_path",
    )
    length(fields[1]) == 128 || error("invalid SHA-512 length in $checksum_path")
    basename(fields[end]) == basename(path) || error(
        "SHA-512 filename does not match source payload: $checksum_path",
    )
    return fields[1]
end

function require_coordinate(dataset, name, expected_size)
    haskey(dataset, name) || error("ECCO source lacks coordinate $name")
    values = Float64.(dataset[name][:])
    size(values) == expected_size || error(
        "ECCO coordinate $name has size $(size(values)), expected $expected_size",
    )
    all(isfinite, values) || error("ECCO coordinate $name contains non-finite values")
    return values
end

expected_start = string(source_month) * "T00:00:00"
expected_end = string(next_month) * "T00:00:00"
reference_longitude = Ref{Any}(nothing)
reference_latitude = Ref{Any}(nothing)
reference_depth = Ref{Any}(nothing)
source_records = NamedTuple[]
variable_sources = Dict{String, String}()

for collection in values(collections)
    source_path = joinpath(raw_directory, collection.filename)
    isfile(source_path) || error("ECCO source payload is missing: $source_path")
    expected_digest = published_sha512(source_path)
    actual_digest = sha512_hex(source_path)
    actual_digest == expected_digest || error(
        "ECCO SHA-512 mismatch for $(collection.filename)",
    )

    NCDataset(source_path, "r") do dataset
        String(dataset.attrib["product_version"]) == "Version 4, Release 4" ||
            error("ECCO source is not V4r4: $(collection.filename)")
        String(dataset.attrib["time_coverage_start"]) == expected_start || error(
            "ECCO source has wrong coverage start: $(collection.filename)",
        )
        String(dataset.attrib["time_coverage_end"]) == expected_end || error(
            "ECCO source has wrong coverage end: $(collection.filename)",
        )
        longitude = require_coordinate(dataset, "longitude", (720,))
        latitude = require_coordinate(dataset, "latitude", (360,))
        all(diff(longitude) .> 0) || error("ECCO longitude is not increasing")
        all(diff(latitude) .> 0) || error("ECCO latitude is not increasing")
        if isnothing(reference_longitude[])
            reference_longitude[] = longitude
            reference_latitude[] = latitude
        else
            longitude == reference_longitude[] || error(
                "ECCO source longitude grids do not match",
            )
            latitude == reference_latitude[] || error(
                "ECCO source latitude grids do not match",
            )
        end

        for variable_name in collection.variables
            haskey(dataset, variable_name) || error(
                "ECCO source $(collection.filename) lacks $variable_name",
            )
            variable = dataset[variable_name]
            expected_shape = variable_name in ("THETA", "SALT", "EVEL", "NVEL") ?
                (720, 360, 50, 1) : (720, 360, 1)
            size(variable) == expected_shape || error(
                "ECCO $variable_name has size $(size(variable)), expected $expected_shape",
            )
            values_array = variable[:]
            finite_values = collect(skipmissing(vec(values_array)))
            isempty(finite_values) && error("ECCO $variable_name contains no valid values")
            all(isfinite, finite_values) || error(
                "ECCO $variable_name contains non-finite valid values",
            )
            variable_sources[variable_name] = source_path
        end

        if any(name -> name in collection.variables, ("THETA", "SALT", "EVEL", "NVEL"))
            depth = require_coordinate(dataset, "Z", (50,))
            all(diff(depth) .< 0) || error("ECCO depth coordinate is not decreasing")
            if isnothing(reference_depth[])
                reference_depth[] = depth
            else
                depth == reference_depth[] || error(
                    "ECCO source depth grids do not match",
                )
            end
        end
    end

    push!(source_records, (
        collection = collection.id,
        filename = collection.filename,
        bytes = filesize(source_path),
        sha512 = actual_digest,
        url = "https://archive.podaac.earthdata.nasa.gov/" *
              "podaac-ops-cumulus-protected/$(collection.id)/$(collection.filename)",
    ))
end

expected_variables = Set(("THETA", "SALT", "EVEL", "NVEL", "SSH", "SIheff", "SIarea"))
Set(keys(variable_sources)) == expected_variables || error(
    "prepared ECCO variable set is incomplete: $(sort!(collect(keys(variable_sources))))",
)

mkpath(output_directory)
for variable_name in sort!(collect(expected_variables))
    source_path = variable_sources[variable_name]
    link_path = joinpath(output_directory, "$(variable_name)_$(filename_label).nc")
    relative_source = relpath(source_path, output_directory)
    if islink(link_path)
        readlink(link_path) == relative_source || error(
            "existing ECCO link points to the wrong source: $link_path",
        )
    elseif ispath(link_path)
        error("refusing to replace existing non-symlink ECCO path: $link_path")
    else
        symlink(relative_source, link_path)
    end
    realpath(link_path) == realpath(source_path) || error(
        "ECCO link does not resolve to its verified source: $link_path",
    )
end

receipt_lines = String[
    "schema=readiesm_ecco_monthly_bundle_v1",
    "dataset=ECCO_V4r4_interpolated_monthly_0p50deg",
    "source_month=$month_label",
    "coverage_start=$expected_start",
    "coverage_end=$expected_end",
    "layout=verified_payloads_with_relative_variable_symlinks",
    "variables=" * join(sort!(collect(expected_variables)), ","),
]
for record in source_records
    prefix = replace(lowercase(record.collection), r"[^a-z0-9]+" => "_")
    append!(receipt_lines, (
        "$(prefix)_filename=$(record.filename)",
        "$(prefix)_bytes=$(record.bytes)",
        "$(prefix)_sha512=$(record.sha512)",
        "$(prefix)_url=$(record.url)",
    ))
end
receipt = join(receipt_lines, "\n") * "\n"
receipt_path = joinpath(output_directory, "ECCO4Monthly_$(filename_label)_RECEIPT.txt")
if isfile(receipt_path)
    read(receipt_path, String) == receipt || error(
        "existing ECCO receipt differs from verified source: $receipt_path",
    )
else
    temporary_path = receipt_path * ".tmp"
    ispath(temporary_path) && error("temporary receipt already exists: $temporary_path")
    open(temporary_path, "w") do io
        write(io, receipt)
    end
    mv(temporary_path, receipt_path)
end

println(
    "ECCO_MONTHLY_BUNDLE_PASS month=$month_label variables=7 sources=4 " *
    "output=$output_directory receipt=$receipt_path",
)
