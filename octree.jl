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
    id_vert::Vector{Int}
    id_child::Union{Vector{Int}, Nothing}
    function Node(id, depth, b_min, b_max)
        ndim = length(b_min)
        new(id, ndim, depth, b_min, b_max, Vector{Int}[], nothing)
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
        N_vert = 0
        depth_init = 1
        vertex = Vector{Vertex}[]
        data = Vector{Float64}[]
        node_root = Node(N_node, depth_init, b_min, b_max)
        new(N_node, N_vert, ndim, depth_init, [node_root], node_root,
            vertex, data, func)
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
            v_lst = bound2vert(node.b_min, node.b_max)
            for v in v_lst
                tree.N_vert += 1
                push!(tree.vertex, v)
                push!(tree.data, tree.func(v))
                push!(node.id_vert, tree.N_vert)
            end
        end
    end
    recursion(tree.node_root)
    println("finish autosplit")
end
function remove_duplicated_vertex!(tree::Tree)
    println("refresing stored data...")
    println("current vertex num is "*string(tree.N_vert))

    # delete vertex and data
    tree.vertex = Vector{Vertex}[]
    tree.data = Vector{Float64}[]
    tree.N_vert = 0
    
    N_vert_1dim = 2^(tree.depth_max-1)+1
    N_vert_whole = N_vert_1dim^tree.ndim
    foot_print = [-1 for i in 1:N_vert_whole] # -1 means unvisited
    b_min_root = tree.node_root.b_min
    b_max_root = tree.node_root.b_max
    size_min = (b_max_root - b_min_root)/2^(tree.depth_max-1)

    visitor_counter = 0

    function recursion(node::Node)
        if node.id_child != nothing
            for id in node.id_child
                recursion(tree.node[id])
            end
        else
            # delete id_vert  
            node.id_vert = Vector{Int}[]
            
            v_lst = bound2vert(node.b_min, node.b_max)
            for v in v_lst
                i = round(Int, (v[1]-b_min_root[1])/size_min[1] + 1)
                j = round(Int, (v[2]-b_min_root[2])/size_min[2] + 1)
                if tree.ndim == 2
                    idx_of_id = (i-1)*N_vert_1dim^1 + (j-1) + 1
                elseif tree.ndim == 3
                    k = round(Int, (v[3]-b_min_root[3])/size_min[3] + 1)
                    idx_of_id = (i-1)*N_vert_1dim^2 + (j-1)*N_vert_1dim + (k-1) + 1
                end

                if foot_print[idx_of_id]==-1 # unvisited
                    visitor_counter += 1
                    foot_print[idx_of_id] = visitor_counter # foot print so that next comer can follow this
                    push!(node.id_vert, visitor_counter)
                else # reach idx where someone else already reached
                    push!(node.id_vert, foot_print[idx_of_id])
                end
                push!(tree.vertex, v)
                push!(tree.data, tree.func(v))
            end

        end
    end
    recursion(tree.node_root)
    tree.N_vert = visitor_counter
    println("reduced vertex num is "*string(tree.N_vert))
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

function show_contour2(tree::Tree; axis_slice=1, pos_slice=0, N_cell=300)
    b_min = tree.node_root.b_min
    b_max = tree.node_root.b_max
    step = (b_max - b_min)./N_cell
    data_plot = zeros(N_cell, N_cell)
    for i in 1:N_cell, j in 1:N_cell
        p_cell = b_min + [i-0.5, j-0.5].*step
        val = evaluate(tree, p_cell)
        data_plot[i, j] = val
    end
    extent = [b_min[1] + 0.5*step[1],
              b_max[1] - 0.5*step[1],
              b_min[2] + 0.5*step[2],
              b_max[2] - 0.5*step[2]]
    PyPlot.imshow(data_plot,
                  extent = extent,
                  interpolation=:bicubic)
    show()

end

function show_contour3(tree::Tree; axis_slice=1, pos_slice=0, N_cell=300)
    # consider right hand 
    if axis_slice == 1
        plot_axis_1 = 2; plot_axis_2 = 3
    elseif axis_slice == 2
        plot_axis_1 = 3; plot_axis_2 = 1
    else
        plot_axis_1 = 1; plot_axis_2 = 2
    end

    b_min_ = tree.node_root.b_min
    b_max_ = tree.node_root.b_max

    b_min = [b_min_[plot_axis_1], b_min_[plot_axis_2]]
    b_max = [b_max_[plot_axis_1], b_max_[plot_axis_2]]
    step = (b_max - b_min)./N_cell

    data_plot = zeros(N_cell, N_cell)
    for i in 1:N_cell, j in 1:N_cell
        p_cell = b_min + [i-0.5, j-0.5].*step
        p_cell_3d = vcat(pos_slice, p_cell)
        val = evaluate(tree, p_cell_3d)
        data_plot[i, j] = val
    end
    extent = [b_min[1] + 0.5*step[1],
              b_max[1] - 0.5*step[1],
              b_min[2] + 0.5*step[2],
              b_max[2] - 0.5*step[2]]
    PyPlot.imshow(data_plot,
                  extent = extent,
                  interpolation=:bicubic)
    show()
end
