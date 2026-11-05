using ReadyESM, Test

@testset "lazy CairoMakie first load" begin
    mktempdir() do directory
        figure = ReadyESM.Figure(size = (240, 160))
        axis = ReadyESM.Axis(ReadyESM._plot_layout(figure, 1, 1))
        ReadyESM.lines!(axis, [0.0, 1.0], [0.0, 1.0])
        path = joinpath(directory, "first_lazy_figure.png")
        ReadyESM.save(path, figure)
        @test isfile(path)
        @test filesize(path) > 0
    end
end
