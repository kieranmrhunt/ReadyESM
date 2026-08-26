"""
Allocate uninitialized storage for a RingGrids `Field` without reading it.

RingGrids 0.3.4 implements `similar(field::F)` by calling the concrete `F`
constructor. That constructor applies `T.(data)` even when `data` already has
element type `T`, so an ordinary `similar` performs an identity broadcast over
the newly allocated, undefined storage. CUDA initcheck correctly reports one
undefined read per element. Constructing through the non-parametric `Field`
inner constructor preserves the exact data and grid types without inspecting
the storage values.
"""
function Base.similar(field::SpeedyWeather.RingGrids.Field)
    return SpeedyWeather.RingGrids.Field(similar(field.data), field.grid)
end
