using LinearAlgebra
using Interpolations
using SpecialFunctions
using PyPlot
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
    b_center = b_min .+ dx[1] .+ dx[2]

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

function needSplitting(node::Node, f)
    ndim = node.ndim
    v_lst = bound2vert(node.b_min, node.b_max)
    f_lst = [f(v) for v in v_lst]

    # TODO fix this by using @eval
    # just for now
    if ndim == 2
        data = zeros(2, 2)
    elseif ndim == 3
        data = zeros(2, 2, 2)
    end

    it = itr(ndim)
    for i = 1:2^ndim
        idc = it()
        if ndim == 2
            data[idc[1]+1, idc[2]+1] = f_lst[i]
        elseif ndim == 3
            data[idc[1]+1, idc[2]+1, idx[3]+1] = f_lst[i]
        end
    end

    # see utils:
    #data = [f_lst[1] f_lst[2]; f_lst[3] f_lst[4]]
    itp = interpolate(data, BSpline(Linear()))

    if ndim == 2
        center_itp = itp(1.5, 1.5) # note: interp starts from 1 
    elseif ndim == 3
        center_itp = itp(1.5, 1.5, 1.5)
    end

    center_real = f(0.5*(node.b_max + node.b_min))

    error = abs(center_itp - center_real)
    #println(error)

    return ~(error<0.02)
end

function auto_split!(tree::Tree, f) # recursive way
    # in the laef, we must add itp
    function recursion(node::Node)
        if needSplitting(node, f)
            split!(tree, node)
            for id in node.id_child
                recursion(tree.node[id])
            end
        else  # if the node is the final decendent, we endow itp to them
            b_min = node.b_min
            b_max = node.b_max
            v_lst = bound2vert(b_min, b_max)
            f_lst = [f(v) for v in v_lst]
            data = [f_lst[1] f_lst[4]; f_lst[2] f_lst[3]]
            itp_ = interpolate(data, BSpline(Linear())) # this raw itp object is useless as it is now

            node.itp = function itp(p)
                p_modif = (p - b_min)./(b_max - b_min) .+ 1
                return itp_(p_modif[1], p_modif[2])
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
            sleep(0.002)
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



sigma = 8
#=
elp(x, y) = x^2/10^2 - y^2/10^2 - 1
f(x) = 0.5*(1 + erf((-elp(x[1], x[2]))/sqrt(2*sigma^2)))
=#
R(a) = [cos(a) -sin(a);
        sin(a) cos(a)]
function sdf(x, a, b)
    x = R(-a)*x
    d = abs.(x) - b
    tmp = [max(d[1], 0.0), max(d[2], 0.0)]
    return norm(tmp) + min(max(d[1], d[2]), 0.0)
end

a = pi/10
b = [60, 15]
f(x) = 0.5*(1 + erf(sdf(x, a, b))/sqrt(2*sigma^2))

t = Tree([-100, -100], [100, 100])
auto_split!(t, f)
show(t)




    

