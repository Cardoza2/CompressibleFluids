module CompressibleFluids

greet() = print("Hello World!")

using Roots
using Plots
using FLOWMath

include("gasproperties.jl")
include("isentropicflow.jl")
include("staticshocks.jl")
include("movingshocks.jl")
include("obliqueshocks.jl")
include("expansionfans.jl")
include("fannoflow.jl")
include("rayleighflow.jl")
include("numericalapproximations.jl")
include("geometry.jl")
include("conversions.jl")


end # module
