module AdaptiveGrid

using LinearAlgebra
using Interpolations

using SpecialFunctions
using StaticArrays

using PyPlot
import Base: show
using JSON

export Tree, Node
export auto_split!, evaluate, pred_standard
export write_json, show, show_contour2, show_contour3
export Linear, Constant # from Interpolations.jl

include("utils.jl")
include("octree.jl")
include("viewer.jl")
include("json_io.jl")

end
