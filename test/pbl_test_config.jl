# Shared by the PBL test modules after their local with(config; kwargs...) helper.
function legacy_pbl_test_config()
    released = R.load_config(joinpath(R.PROJECT_ROOT, "config/production.yml"); check_input_files=false)
    @test released.atmosphere_vertical_diffusion == :holtslag_boville_k_profile_v1
    @test released.atmosphere_pbl_diffusivity_scale == 1.0
    @test released.atmosphere_pbl_mix_momentum
    @test !released.atmosphere_pbl_mix_heat && !released.atmosphere_pbl_mix_humidity
    # These scalar and refusal tests exercise the old zero-mixing starting point.
    # Their fixture must remain independent of changes to the released defaults.
    return with(released; atmosphere_vertical_diffusion=:speedyweather_0_21_1_zero_operator_control,
        atmosphere_pbl_diffusivity_scale=1.0, atmosphere_pbl_mix_momentum=false,
        atmosphere_pbl_mix_heat=true, atmosphere_pbl_mix_humidity=true)
end

function pbl_fixture_yaml()
    # Remove keys before adding explicit test settings, so YAML has no duplicates.
    keys = ("atmosphere_pbl_diffusivity_scale:", "atmosphere_pbl_mix_momentum:",
            "atmosphere_pbl_mix_heat:", "atmosphere_pbl_mix_humidity:")
    lines = filter(line -> !any(key -> startswith(line, key), keys),
                   readlines(joinpath(R.PROJECT_ROOT, "config/production.yml")))
    return lines
end
