include("octree.jl")
using JSON

## json writer 
function write_json(tree::Tree, filename)
    dict = Dict("N_node"=>tree.N_node,
             "N_vert"=>tree.N_vert,
             "ndim"=>tree.ndim,
             "depth_max"=>tree.depth_max,
             "node"=>Dict(string(i)=>convert_to_dict(tree.node[i]) for i in 1:length(tree.node)),
             "vertex"=>Dict(string(i)=>tree.vertex[i] for i in 1:length(tree.vertex)),
             "data"=>Dict(string(i)=>tree.data[i] for i in 1:length(tree.data)))
    j = JSON.json(dict)
    open(filename, "w") do f
        write(f, j)
    end
end

function convert_to_dict(node::Node)
    d = Dict([("id", node.id),
              ("ndim", node.ndim),
              ("depth", node.depth),
              ("b_min", node.b_min),
              ("b_max", node.b_max),
              ("id_vert", node.id_vert),
              ("id_child", node.id_child)])
    return d
end

## json reader (constructor)
function Tree(filename)
    text = read(filename, String)
    dict = JSON.parse(text)
    N_node = dict["N_node"]
    N_vert = dict["N_vert"]
    ndim = dict["ndim"]
    depth_max = dict["depth_max"]
    node = [Node(dict["node"][string(i)]) for i in 1:N_node]
    node_root = node[1]
    vertex = [convert(SVector{ndim, Float64}, dict["vertex"][string(i)]) for i in 1:N_vert]
    data = [dict["data"][string(i)] for i in 1:N_vert]

    # check
    length(node)!=N_node && error("???")
    length(vertex)!=N_vert && error("???")
    return Tree(N_node, N_vert, ndim, depth_max, node, node_root, vertex, data)
end

function Node(dict)
    id = dict["id"]
    ndim = dict["ndim"]
    depth = dict["depth"]

    if ndim == 2
        b_min = SVector(dict["b_min"][1], dict["b_min"][2])
        b_max = SVector(dict["b_max"][1], dict["b_max"][2])
    else
        b_min = SVector(dict["b_min"][1], dict["b_min"][2], dict["b_min"][3])
        b_max = SVector(dict["b_max"][1], dict["b_max"][2], dict["b_max"][3])
    end
    
    id_vert = dict["id_vert"]
    if id_vert != nothing
        id_vert = convert(Vector{Int}, id_vert)
    end

    id_child = dict["id_child"]
    if id_child != nothing
        id_child = convert(Vector{Int}, id_child)
    end

    return Node(id, ndim, depth, b_min, b_max, id_vert, id_child)
end

