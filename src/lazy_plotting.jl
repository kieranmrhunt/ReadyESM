const _CAIRO_MAKIE_PACKAGE_ID = Base.PkgId(
    Base.UUID("13f3f980-e62b-5c42-98c6-ff1f3baf88f0"),
    "CairoMakie",
)
const _CAIRO_MAKIE_MODULE = Ref{Union{Nothing, Module}}(nothing)

"""
Load CairoMakie only when a plotting entry point is actually called.

Production integrations save NetCDF and render figures in separate processes,
so eagerly importing the complete plotting stack makes every model build pay an
otherwise unused startup and memory cost. `Base.require` returns the package
module without injecting its exports into ReadyESM; the narrow wrappers below
preserve the existing internal plotting calls.
"""
function _cairo_makie()
    module_ = _CAIRO_MAKIE_MODULE[]
    if isnothing(module_)
        module_ = Base.require(_CAIRO_MAKIE_PACKAGE_ID)
        _CAIRO_MAKIE_MODULE[] = module_
    end
    return module_
end

function _cairo_makie_call(name::Symbol, args...; kwargs...)
    function_ = getproperty(_cairo_makie(), name)
    return Base.invokelatest(function_, args...; kwargs...)
end

for name in (
    :Figure,
    :Axis,
    :Colorbar,
    :axislegend,
    :lines!,
    :scatter!,
    :heatmap!,
    :save,
)
    @eval begin
        function $(name)(args...; kwargs...)
            return _cairo_makie_call($(QuoteNode(name)), args...; kwargs...)
        end
    end
end

# Figure layout indexing dispatches through a `Base.getindex` method that Makie
# defines when it loads. A plotting function compiled before that first load
# otherwise sees the older world and fails even though constructors above use
# `invokelatest`. Keep layout lookup behind the same explicit world-age boundary.
_plot_layout(figure, indices...) = Base.invokelatest(getindex, figure, indices...)
