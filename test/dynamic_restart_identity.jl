@testset "dynamic restart identity" begin
    @test ReadyESM._require_rrtmgp_surface_albedo_restart((;
        direct_surface_albedo = [0.1], diffuse_surface_albedo = [0.2])) === nothing
    @test_throws r"predates RRTMGP surface-albedo restart state" begin
        ReadyESM._require_rrtmgp_surface_albedo_restart((;
            direct_surface_albedo = [0.1]))
    end
    @test_throws r"predates RRTMGP surface-albedo restart state" begin
        ReadyESM._require_rrtmgp_surface_albedo_restart((;))
    end
    config = ReadyESM.ExperimentConfig()
    changed_config = ReadyESM.ExperimentConfig(name = "identity-test")
    @test ReadyESM._dynamic_config_identity_hash(config) ==
          ReadyESM._dynamic_config_identity_hash(config)
    @test ReadyESM._dynamic_config_identity_hash(config) ==
          ReadyESM._dynamic_config_identity_hash(changed_config)
    for field in ReadyESM._DYNAMIC_RESTART_RUN_CONTROL_FIELDS
        changed = ReadyESM.ExperimentConfig(; merge(
            (; (name => getfield(config, name) for name in fieldnames(typeof(config)))...),
            (; field => field in (:duration_days, :analysis_start_day) ? 20.0 :
                field in (:scored_outputs, :sea_ice_branch_observer) ? true : "changed"),
        )...)
        @test ReadyESM._dynamic_config_identity_hash(config) ==
              ReadyESM._dynamic_config_identity_hash(changed)
    end
    @test ReadyESM._dynamic_config_identity_hash(config) !=
          ReadyESM._dynamic_config_identity_hash(ReadyESM.ExperimentConfig(
              atmosphere_divergence_hyperdiffusion_hours = 4.0,
          ))

    # Absolute and relative names of an input inside this checkout have the
    # same logical identity; an external path remains distinct.
    project = ReadyESM.PROJECT_ROOT
    absolute = ReadyESM.ExperimentConfig(
        ecco_initial_conditions_directory = joinpath(project, "data", "ecco_fixture"),
        forcing = ReadyESM.ForcingConfig(
            ozone_path = joinpath(project, "data", "ozone_fixture.nc"),
        ),
    )
    relative = ReadyESM.ExperimentConfig(
        ecco_initial_conditions_directory = "data/ecco_fixture",
        forcing = ReadyESM.ForcingConfig(ozone_path = "data/ozone_fixture.nc"),
    )
    external = ReadyESM.ExperimentConfig(
        ecco_initial_conditions_directory = joinpath(dirname(project), "outside", "ecco_fixture"),
        forcing = absolute.forcing,
    )
    @test ReadyESM._dynamic_restart_portable_path(absolute.ecco_initial_conditions_directory) ==
          "project:data/ecco_fixture"
    @test ReadyESM._dynamic_config_identity_hash(absolute) ==
          ReadyESM._dynamic_config_identity_hash(relative)
    @test ReadyESM._dynamic_config_identity_hash(absolute) !=
          ReadyESM._dynamic_config_identity_hash(external)
    @test ReadyESM._dynamic_config_identity_hash(absolute) !=
          ReadyESM._dynamic_config_identity_hash(ReadyESM.ExperimentConfig(
              ecco_initial_conditions_directory = absolute.ecco_initial_conditions_directory,
              forcing = ReadyESM.ForcingConfig(ozone_path = "data/other_ozone.nc"),
          ))

    # `repr(DateTime)` changes with the caller's Dates import. The identity
    # must be stable across fresh Julia processes with both import contexts.
    julia = Base.julia_cmd()
    project_arg = ReadyESM.PROJECT_ROOT
    bare = readchomp(`$julia --startup-file=no --project=$project_arg -e "using ReadyESM; print(ReadyESM._dynamic_config_identity_hash(ReadyESM.ExperimentConfig()))"`)
    with_dates = readchomp(`$julia --startup-file=no --project=$project_arg -e "using ReadyESM, Dates; print(ReadyESM._dynamic_config_identity_hash(ReadyESM.ExperimentConfig()))"`)
    @test bare == with_dates == ReadyESM._dynamic_config_identity_hash(config)

    model = Ref{Nothing}(nothing)
    fixture_commit = repeat("a", 40)
    withenv("READYESM_SOURCE_COMMIT" => fixture_commit) do
        ReadyESM._bind_dynamic_restart_identity!(model, config)
        identity = ReadyESM._dynamic_restart_identity(model)
        @test identity.schema == "readiesm-dynamic-restart-identity-v1"
        @test identity.config_sha256 == ReadyESM._dynamic_config_identity_hash(config)
        @test identity.source_commit == fixture_commit
        @test occursin(r"^[0-9a-f]{64}$", identity.manifest_sha256)
        @test ReadyESM._validate_dynamic_restart_identity!(
            model,
            (; restart_identity = identity),
        )

        changed_config_hash = merge(identity, (; config_sha256 = repeat("0", 64)))
        changed_commit = merge(identity, (; source_commit = repeat("b", 40)))
        changed_manifest = merge(identity, (; manifest_sha256 = repeat("0", 64)))
        for changed_identity in (changed_config_hash, changed_commit, changed_manifest)
            @test_throws ErrorException ReadyESM._validate_dynamic_restart_identity!(
                model,
                (; restart_identity = changed_identity),
            )
        end
        @test_throws ErrorException ReadyESM._validate_dynamic_restart_identity!(
            model,
            (;),
        )
        @test ReadyESM.with_dynamic_restart_identity_override() do
            !ReadyESM._validate_dynamic_restart_identity!(
                model,
                (; restart_identity = changed_manifest),
            )
        end
    end
end
