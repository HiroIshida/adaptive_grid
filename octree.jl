using LinearAlgebra
import Interpolations
using SpecialFunctions
using PyPlot
import Base: show
using NearestNeighbors
include("utils.jl")

const Vertex = Vector{Float64}

mutable struct Node
    id::Int
    ndim::Int
    depth::Int
    b_min
    b_max
    id_vert::Union{Vector{Int}, Nothing}
    id_child::Union{Vector{Int}, Nothing}
    function Node(id, depth, b_min, b_max)
        ndim = length(b_min)
        new(id, ndim, depth, b_min, b_max, nothing, nothing)
    end
end

function show(node::Node; color=:r)
    v_lst = bound2vert(node.b_min, node.b_max)
    if node.ndim == 2
        lst_pair = [[1, 2], [2, 4], [4, 3], [3, 1]]
        for pair in lst_pair 
            idx1 = pair[1]; idx2 = pair[2]
            x = [v_lst[idx1][1], v_lst[idx2][1]]
            y = [v_lst[idx1][2], v_lst[idx2][2]]
            PyPlot.plot(x, y, color)
        end
    elseif node.ndim == 3
        lst_pair = [[1, 3], [3, 4], [4, 2], [2, 1], 
                    [1, 5], [2, 6], [4, 8], [3, 7],
                    [5, 7], [7, 8], [8, 6], [6, 5]]
        for pair in lst_pair 
            idx1 = pair[1]; idx2 = pair[2]
            x = [v_lst[idx1][1], v_lst[idx2][1]]
            y = [v_lst[idx1][2], v_lst[idx2][2]]
            z = [v_lst[idx1][3], v_lst[idx2][3]]
            PyPlot.plot3D(x, y, z, color)
        end
    else
        error("is not supported")
    end
end

mutable struct Tree
    # member variable
    N_node::Int
    N_vert::Int
    ndim::Int
    depth_max::Int
    node::Vector{Node}
    node_root::Node
    vertex::Vector{Vertex}
    data::Vector{Float64} # data stored in vertex
    func#function

    function Tree(b_min, b_max, func)
        ndim = length(b_min)
        N_node = 1
        N_vert = 2^ndim
        depth_init = 1
        v_lst = bound2vert(b_min, b_max)
        f_lst = [func(v) for v in v_lst]
        node_root = Node(N_node, depth_init, b_min, b_max)
        new(N_node, N_vert, ndim, depth_init, [node_root], node_root, v_lst, f_lst, func)
    end
end

function split!(tree::Tree, node::Node)
    # id of new node 
    id_child = [i for i in 1:2^(tree.ndim)] .+ tree.N_node
    node.id_child = id_child

    # generate new nodes
    b_min = node.b_min
    b_max = node.b_max
    dif = b_max - b_min
    dx = bound2dx(b_min, b_max)*0.5
    b_center = (b_min + b_max)*0.5

    for i in 0:2^tree.ndim-1
        add = [0.0 for n in 1:tree.ndim]
        for dim in 1:tree.ndim
            if mod(div(i, 2^(dim-1)), 2) == 1
                add += dx[dim]
            end
        end

        # edit new node and corresponding vertex in tree
        id_new = tree.N_node+1
        b_min_new = b_min+add
        b_max_new = b_center+add # not b_max + add
        depth_new = node.depth+1
        node_new = Node(tree.N_node+1, depth_new, b_min_new, b_max_new)
        push!(tree.node, node_new)
        tree.N_node += 1
        tree.depth_max < depth_new && (tree.depth_max = depth_new)
    end
end


function auto_split!(tree::Tree, predicate) 
    # recusive split based on the boolean returned by predicate
    # perdicate: Node →  bool
    # interpolation objecet is endowed with each terminal nodes hh
    function recursion(node::Node)
        if predicate(node)
            split!(tree, node)
            for id in node.id_child
                recursion(tree.node[id])
            end
        else
            error("ここでvertexを挿入しろ")
        end
    end
    recursion(tree.node_root)
    println("finish autosplit")
end

function remove_duplicated_vertex!(tree::Tree; map_cache=nothing, ids_cache=nothing)
    println("start vertex reductoin")
    # first re-label the indices.
    # for example if S1 = [1, 4, 6], S2 = [2, 3, 7], S3 =[5, 8] are duplicated
    # label them i1=1 i2=2 i3=3. Then, make a map from S -> i
    # potentially dangerous operation

    if (map_cache==nothing) | (ids_cache==nothing) 
        map, valid_ids = vertex_reduction(tree.vertex)
    else # if use cache
        map = map_cache
        valid_ids = ids_cache
    end

    vertex_new = Vertex[]
    data_new = Float64[]
    for id in valid_ids
        push!(vertex_new, tree.vertex[id])
        push!(data_new, tree.data[id])
    end
    tree.N_vert = length(vertex_new)
    tree.vertex = vertex_new
    tree.data = data_new

    function recursion(node::Node)
        if node.id_child!=nothing
            for id in node.id_child
                recursion(tree.node[id])
            end
        else
            node.id_vert = [map[i] for i in node.id_vert]
        end
    end
    recursion(tree.node_root)
    println("end vertex reduction")
    return map, valid_ids
end

function show(tree::Tree)
    function recursion(node::Node)
        if node.id_child!=nothing
            for id in node.id_child
                recursion(tree.node[id])
            end
        else
            show(node; color=:r)
        end
    end
    recursion(tree.node_root)
    println("finish show")
end

function evaluate(tree::Tree, q)
    node = tree.node_root
    while(true)
        if node.id_child == nothing
            data_itp = form_data_cubic(tree.data[node.id_vert], tree.ndim)
            itp = interpolate(data_itp, BSpline(Linear())) #raw
            q_modif = (q - node.b_min)./(node.b_max - node.b_min) .+ 1
            if tree.ndim == 2
                return itp(q_modif[1], q_modif[2])
            elseif tree.ndim == 3
                return itp(q_modif[1], q_modif[2], q_modif[3])
            else
                error("not supported")
            end
        end
        idx = whereami(q, node.b_min, node.b_max)
        id_next = node.id_child[idx]
        node = tree.node[id_next]
    end
end

function evaluate_(tree::Tree, q)
    node = tree.node_root
    while(true)
        node.id_child == nothing && return node.itp(q)
        idx = whereami(q, node.b_min, node.b_max)
        id_next = node.id_child[idx]
        node = tree.node[id_next]
    end
end


function pred_standard(node::Node, f, ε, n_grid, itp_method)
    # choose interpolation method: itp_method
    # curretnly available methods are "Linear()" and "Constant()"
    # apparently, if you choose "Constant", much more cells are required.
    
    ndim = node.ndim
    v_lst = bound2vert(node.b_min, node.b_max)
    f_lst = [f(v) for v in v_lst]
    data = form_data_cubic(f_lst, ndim)
    itp = interpolate(data, BSpline(itp_method))

    # evaluate the interpolation error for many points in the cell 
    # and only care about maximum error 
    points_eval = grid_points(n_grid, node.b_min, node.b_max)
    max_error = -Inf
    for p in points_eval
        p_reg = (p - node.b_min)./(node.b_max - node.b_min) .+ 1
        val_itp = (ndim == 2 ? itp(p_reg[1], p_reg[2]) : itp(p_reg[1], p_reg[2], p_reg[3]))
        val_real = f(p)
        error = abs(val_real - val_itp)
        if max_error < error
            max_error = error
        end
    end
    return ~(max_error<ε)
end

