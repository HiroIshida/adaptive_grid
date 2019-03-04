push!(LOAD_PATH, "./NearestNeighbors.jl")
using NearestNeighbors
include("utils.jl")
data = grid_points(10)
list = grid_points(10)
for v in list
    push!(data, v)
end

m = zeros(3, length(data))
for i in 1:(length(data))
    m[:, i] = [data[i][1], data[i][2], data[i][3]]
end

kdtree = KDTree(m; leafsize=10)
inrange(kdtree, list[100], 0.0001, true)

