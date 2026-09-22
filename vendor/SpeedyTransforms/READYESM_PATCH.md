# ReadyESM patches

This copy starts from SpeedyTransforms 0.1.5. The upstream stable-pointer graph
cache change (SpeedyWeather commit `54d500f4d5eb006cc83998978509568d7efa4818`)
avoids a new cache entry whenever a CuArray view creates a new Julia wrapper.

ReadyESM additionally includes both Fourier scratch arrays, their dimensions
and strides in that key. The public transform API accepts explicit scratch
storage. Replaying a graph captured for other scratch buffers otherwise writes
to or reads from those previous buffers, silently returning the wrong result.
The executable now retains its input/output arrays, preventing allocation
reuse while a cached graph still contains their pointers. Cache clearing waits
for graph execution to complete before releasing those owners.

The forward and inverse scratch-switch regression fails on the previous
field-pointer-only implementation. See `test/cuda_graph_cache.jl` for the GPU
regression, including view reuse and buffer lifetime checks. This repair does
not by itself establish the cause of the separate coupled CUDA restart fault.
