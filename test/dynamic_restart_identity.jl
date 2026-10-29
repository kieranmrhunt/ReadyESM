@testset "dynamic restart identity" begin
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
                field == :scored_outputs ? true : "changed"),
        )...)
        @test ReadyESM._dynamic_config_identity_hash(config) ==
              ReadyESM._dynamic_config_identity_hash(changed)
    end
    @test ReadyESM._dynamic_config_identity_hash(config) !=
          ReadyESM._dynamic_config_identity_hash(ReadyESM.ExperimentConfig(
              atmosphere_divergence_hyperdiffusion_hours = 4.0,
          ))

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
