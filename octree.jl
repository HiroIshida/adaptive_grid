using LinearAlgebra
using Interpolations
using SpecialFunctions
using PyPlot
using Test
import Base: show
include("utils.jl")

mutable struct Node
    id::Int
    ndim::Int
    b_min
    b_max
    id_child::Union{Vector{Int}, Nothing}
    itp # type?
    function Node(id, b_min, b_max)
        ndim = length(b_min)
        new(id, ndim, b_min, b_max, nothing, nothing)
    end
end

function show(node::Node; color=:r)
    v_lst_ = bound2vert(node.b_min, node.b_max)
    v_lst = [v_lst_[1], v_lst_[2], v_lst_[4], v_lst_[3]]
    for n = 1:4
        if n!=4
            x = [v_lst[n][1], v_lst[n+1][1]]
            y = [v_lst[n][2], v_lst[n+1][2]]
        else
            x = [v_lst[4][1], v_lst[1][1]]
            y = [v_lst[4][2], v_lst[1][2]]
        end
        PyPlot.plot(x, y, color)
    end
end

mutable struct Tree
    # member variable
    N::Int
    ndim::Int
    node::Vector{Node}
    node_root::Node

    function Tree(b_min, b_max)
        id_node = 1
        ndim = length(b_min)
        node_root = Node(id_node, b_min, b_max)
        new(id_node, ndim, [node_root], node_root)
    end
end

function split!(tree::Tree, node::Node)
    # edit node
    leaves = [i for i in 1:2^(tree.ndim)] .+ tree.N
    node.id_child = leaves

    # edit tree
    b_min = node.b_min
    b_max = node.b_max
    dif = b_max - b_min

    dx = bound2dx(b_min, b_max)*0.5
    b_center = (b_min + b_max)*0.5

    # seems complicated but...
    for i in 0:2^tree.ndim-1
        add = [0.0 for n in 1:tree.ndim]
        for dim in 1:tree.ndim
            if mod(div(i, 2^(dim-1)), 2) == 1
                add += dx[dim]
            end
        end
        node_new = Node(tree.N+1, b_min+add, b_center+add)
        push!(tree.node, node_new)
        tree.N += 1
    end
end

function pred_simplest(node::Node, f, ε)
    ndim = node.ndim
    v_lst = bound2vert(node.b_min, node.b_max)
    f_lst = [f(v) for v in v_lst]
    data = form_data_cubic(f_lst, ndim)

    itp = interpolate(data, BSpline(Linear()))

    # TODO dirty dirty dirty dirty
    if ndim == 2
        center_itp = itp(1.5, 1.5) # note: interp starts from 1 
    elseif ndim == 3
        center_itp = itp(1.5, 1.5, 1.5)
    end

    center_real = f((node.b_max + node.b_min)/2)
    error = abs(center_itp - center_real)

    return ~(error<ε)
end

function auto_split!(tree::Tree, f, predicate) 
    # recusive split based on the boolean returned by predicate
    # perdicate: Node →  bool
    # interpolation objecet is endowed with each terminal nodes hh
    
    function recursion(node::Node)
        if predicate(node)
            split!(tree, node)
            for id in node.id_child
                recursion(tree.node[id])
            end
        else  # if terminal node
            b_min = node.b_min
            b_max = node.b_max
            v_lst = bound2vert(b_min, b_max)
            f_lst = [f(v) for v in v_lst]
            data = form_data_cubic(f_lst, tree.ndim)
            itp_ = interpolate(data, BSpline(Linear())) # this raw itp object is useless as it is now


            node.itp = function itp(p)
                p_modif = (p - b_min)./(b_max - b_min) .+ 1
                # TODO dirty, dirty, dirty!!
                if tree.ndim == 2
                    return itp_(p_modif[1], p_modif[2])
                else tree.ndim == 3
                    return itp_(p_modif[1], p_modif[2], p_modif[3])
                end
            end
        end
    end
    recursion(tree.node_root)
    println("finish autosplit")
end

function show(tree::Tree)
    function recursion(node::Node)
        if node.id_child!=nothing
            for id in node.id_child
                recursion(tree.node[id])
            end
        else
            show(node; color=:r)
            #sleep(0.002)
        end
    end
    recursion(tree.node_root)
    println("finish show")
end

function search_idx(tree::Tree, q)
    node = tree.node_root
    while(true)
        if node.id_child == nothing
            return node.id
        end
        idx = whereami(q, node.b_min, node.b_max)
        id_next = node.id_child[idx]
        node = tree.node[id_next]
    end
end

function evaluate(tree::Tree, q)
    id = search_idx(tree, q)
    node = tree.node[id]
    return node.itp(q)
end








    

