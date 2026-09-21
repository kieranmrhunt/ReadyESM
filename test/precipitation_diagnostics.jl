@testset "Native-precision precipitation component checks" begin
    closes = ReadyESM._precipitation_components_close
    # Observed production column: the old absolute-only check rejects this
    # 0.31-native-epsilon relative discrepancy after kg m⁻² s⁻¹ conversion.
    total = [0.0030452770261035766]
    convection = [0.0010029024224422756]
    large_scale = [0.0020423744899744634]
    @test abs(only(total - convection - large_scale)) > 1e-10
    @test closes(total, convection + large_scale)
    @test closes([0.0], [0.0])
    @test closes([5e-11], [0.0])  # Retain the existing near-zero floor.
    @test !closes([2e-10], [0.0])
    @test !closes(total .* (1 + 1e-6), convection + large_scale)
    @test !closes(total, convection)  # Missing large-scale rain.
    @test !closes(total, 2convection + large_scale)
    @test !closes([NaN], [NaN])
    @test !closes([Inf], [Inf])
    @test !closes([1.0], [1.0, 1.0])
    @test !closes([1.0, 1.0], [1.0])
    for rate in Float32[1e-6, 1e-5, 1e-4, 1e-3]
        a, b = rate, rate / 3
        rounded_total = [1000Float64(a + b)]
        component_sum = [1000(Float64(a) + Float64(b))]
        @test closes(rounded_total, component_sum)
        @test !closes(rounded_total .* 1.001, component_sum)
    end
end
