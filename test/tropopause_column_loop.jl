using ReadyESM, Test

@testset "PA2 tropopause is diagnosed for every column" begin
    # Bottom-to-top levels. The first two columns qualify at different layers;
    # the third has no qualifying layer and legitimately falls back to the top.
    pressure = repeat(reshape([90_000.0, 70_000.0, 45_000.0, 30_000.0, 10_000.0], 5, 1), 1, 3)
    height = repeat(reshape([0.0, 2_000.0, 5_000.0, 8_000.0, 12_000.0], 5, 1), 1, 3)
    temperature = hcat(
        [300.0, 280.0, 270.0, 269.0, 268.0],
        [300.0, 280.0, 270.0, 250.0, 249.0],
        [300.0, 280.0, 270.0, 250.0, 220.0],
    )
    @test ReadyESM._tropopause_layers(temperature, pressure, height) == [3, 4, 5]
    @test ReadyESM._tropopause_layers(temperature[:, [2, 1, 3]], pressure, height) == [4, 3, 5]
end

println("PA2_TROPOPAUSE_COLUMN_LOOP_PASS")
