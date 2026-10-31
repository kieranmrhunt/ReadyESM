const CLIMATE_WRITER_TEST_TMP = mktempdir()
ENV["CLIMATE_PROFILE_TEST_WORK"] = joinpath(CLIMATE_WRITER_TEST_TMP, "profile-output")
ENV["CLIMATE_ENERGY_TEST_WORK"] = joinpath(CLIMATE_WRITER_TEST_TMP, "energy-output")
atexit(() -> rm(CLIMATE_WRITER_TEST_TMP; recursive = true, force = true))

module MonthlyWriterCoreTests
include(joinpath(@__DIR__, "monthly_writer_core.jl"))
end

module MonthlyWriterGeometryTests
include(joinpath(@__DIR__, "monthly_writer_geometry.jl"))
end

module MonthlyWriterIOTests
include(joinpath(@__DIR__, "monthly_writer_io.jl"))
end

module MonthlyWriterBindingTests
include(joinpath(@__DIR__, "monthly_writer_binding.jl"))
end

module MonthlyWriterCustomDimsTests
include(joinpath(@__DIR__, "monthly_writer_custom_dims.jl"))
end

module MonthlyWriterSamplingTests
include(joinpath(@__DIR__, "monthly_writer_sampling.jl"))
end

module ExtendedFieldsTests
include(joinpath(@__DIR__, "test_fields.jl"))
end

module ExtendedFieldsSamplingTests
include(joinpath(@__DIR__, "test_sample.jl"))
end

module ExtendedFieldsIOTests
include(joinpath(@__DIR__, "fields_io.jl"))
end

module ClimateProfilesTests
include(joinpath(@__DIR__, "test_profiles.jl"))
end

module ClimateProfilesIOTests
include(joinpath(@__DIR__, "profiles_io.jl"))
end

module ClimateEnergyTests
include(joinpath(@__DIR__, "test_energy.jl"))
end

module ClimateEnergyIOTests
include(joinpath(@__DIR__, "energy_io.jl"))
end

module CMIPFieldsTests
include(joinpath(@__DIR__, "test_cmip_fields.jl"))
end

module CMIPOceanFieldsTests
include(joinpath(@__DIR__, "test_ocean_fields.jl"))
end

module CMIPOverturningTests
include(joinpath(@__DIR__, "test_overturning.jl"))
end

module CMIPOverturningGridTests
include(joinpath(@__DIR__, "test_sample_overturning_grids.jl"))
end
