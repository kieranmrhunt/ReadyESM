@testset "precompile workload mode" begin
    for (raw, mode) in (
        "" => :bounded, "bounded" => :bounded, "default" => :bounded,
        "1" => :bounded, "true" => :bounded, "yes" => :bounded, "on" => :bounded,
        "0" => :off, "false" => :off, "no" => :off, "off" => :off, "none" => :off,
        "full" => :full, "all" => :full, "coupled" => :full,
    )
        @test ReadyESM._precompile_mode(raw) === mode
        @test ReadyESM._precompile_mode("  " * uppercase(raw) * "  ") === mode
    end
    for typo in ("of", "ful", "disable", "unexpected")
        @test_throws ArgumentError ReadyESM._precompile_mode(typo)
    end
    withenv("READYESM_PRECOMPILE_WORKLOAD_MODE" => "off",
            "READYESM_PRECOMPILE_WORKLOAD" => "full") do
        @test ReadyESM._precompile_mode() === :off
    end
    withenv("READYESM_PRECOMPILE_WORKLOAD_MODE" => nothing,
            "READYESM_PRECOMPILE_WORKLOAD" => "full") do
        @test ReadyESM._precompile_mode() === :full
    end
    withenv("READYESM_PRECOMPILE_WORKLOAD_MODE" => "ful",
            "READYESM_PRECOMPILE_WORKLOAD" => "off") do
        @test_throws ArgumentError ReadyESM._precompile_mode()
    end
end
