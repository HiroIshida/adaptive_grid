using StaticArrays
using JSON


mutable struct Node{N}
    id::Int
    ndim::Int
    depth::Int
    b_min::SVector{N, Float64}
    b_max::SVector{N, Float64}
    id_vert::Vector{Int}
    id_child::Union{Vector{Int}, Nothing}
    function Node(id, depth, b_min, b_max)
        ndim = length(b_min)
        new{ndim}(id, ndim, depth, b_min, b_max, Vector{Int}[], nothing)
    end
end

function convert_to_dict(node::Node)
    d = Dict([("id", node.id),
              ("ndim", node.ndim),
              ("b_min", node.b_min),
              ("b_max", node.b_max),
              ("id_child", node.id_child)])
    return d
end

n = Node(1, 1, [0, 0], [1, 1])
dict = Dict(string(i) => convert_to_dict(n) for i in 1:100000)

#j = JSON.json(n)
a="{\"menu\": {
         \"id\": \"file\",
         \"value\": \"File\",
         \"popup\": {
           \"menuitem\": [
             {\"value\": \"New\", \"onclick\": \"CreateNewDoc()\"},
             {\"value\": \"Open\", \"onclick\": \"OpenDoc()\"},
             {\"value\": \"Close\", \"onclick\": \"CloseDoc()\"}
           ]
         }
       }}
       "

