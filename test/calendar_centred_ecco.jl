@testset "Calendar-centred ECCO launch bundle" begin
    include(joinpath(
        ReadyESM.PROJECT_ROOT,
        "scripts",
        "verify_ecco_launch_bundle.jl",
    ))

    @test config.ocean_initial_conditions == :ecco_monthly_full_state
    @test (config.ecco_year, config.ecco_month) == (1993, 1)
    @test config.ecco_initial_conditions_directory == joinpath(
        ReadyESM.PROJECT_ROOT,
        "data",
        "ecco_calendar_centred_1993_01",
    )

    mktempdir() do directory
        qualified = joinpath(directory, "qualified")
        other = joinpath(directory, "other")
        mkpath(qualified)
        mkpath(other)
        @test requires_calendar_centred_ecco_verification(
            "production",
            other,
            qualified,
        )
        @test requires_calendar_centred_ecco_verification(
            "renamed-production-copy",
            qualified,
            qualified,
        )
        @test !requires_calendar_centred_ecco_verification(
            "experiment",
            other,
            qualified,
        )
    end

    function write_test_bundle(directory, manifest_path)
        mkpath(directory)
        entries = Pair{String, String}[]
        for filename in sort!(collect(CALENDAR_CENTRED_ECCO_FILES))
            path = joinpath(directory, filename)
            write(path, "calendar-centred test payload for $filename\n")
            push!(entries, filename => _calendar_centred_ecco_sha256(path))
        end
        open(manifest_path, "w") do stream
            for (filename, digest) in entries
                println(stream, "$digest  $filename")
            end
        end
        return entries
    end

    mktempdir() do directory
        bundle = joinpath(directory, "bundle")
        manifest = joinpath(directory, "manifest.sha256")
        entries = write_test_bundle(bundle, manifest)
        @test isnothing(redirect_stdout(devnull) do
            verify_calendar_centred_ecco_bundle(
                bundle,
                manifest;
                require_caches = true,
            )
        end)

        first_filename = first(entries).first
        open(joinpath(bundle, first_filename), "a") do stream
            write(stream, "corruption\n")
        end
        @test_throws ErrorException verify_calendar_centred_ecco_bundle(
            bundle,
            manifest,
        )
    end

    mktempdir() do directory
        bundle = joinpath(directory, "bundle")
        manifest = joinpath(directory, "manifest.sha256")
        write_test_bundle(bundle, manifest)
        for filename in CALENDAR_CENTRED_ECCO_CACHE_FILES
            rm(joinpath(bundle, filename))
        end
        @test isnothing(redirect_stdout(devnull) do
            verify_calendar_centred_ecco_bundle(bundle, manifest)
        end)
        @test_throws ErrorException verify_calendar_centred_ecco_bundle(
            bundle,
            manifest;
            require_caches = true,
        )
        filename = first(CALENDAR_CENTRED_ECCO_CACHE_FILES)
        write(joinpath(bundle, filename), "partial cache\n")
        @test_throws ErrorException verify_calendar_centred_ecco_bundle(
            bundle,
            manifest,
        )
    end

    mktempdir() do directory
        bundle = joinpath(directory, "bundle")
        manifest = joinpath(directory, "manifest.sha256")
        write_test_bundle(bundle, manifest)
        write(joinpath(bundle, "unexpected.nc"), "unexpected\n")
        @test_throws ErrorException verify_calendar_centred_ecco_bundle(
            bundle,
            manifest,
        )
    end

    mktempdir() do directory
        bundle = joinpath(directory, "bundle")
        manifest = joinpath(directory, "manifest.sha256")
        write_test_bundle(bundle, manifest)
        write(manifest, "malformed\n")
        @test_throws ErrorException verify_calendar_centred_ecco_bundle(
            bundle,
            manifest,
        )
    end
end
