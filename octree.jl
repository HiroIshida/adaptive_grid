using LinearAlgebra
using Interpolations
using PyPlot

function bound2vert(bmin, bmax)
    dif = bmax - bmin
    dx = [dif[1], 0]
    dy = [0, dif[2]]

    v1 = bmin
    v2 = bmin + dx
    v3 = bmin + dx + dy
    v4 = bmin + dy
    v_lst = [v1, v2, v3, v4]
    return v_lst
end


mutable struct Node
    id_node::Int
    b_min
    b_max
    hasLeaf::Bool
    leaf::Union{Vector{Int}, Nothing}
    function Node(id_node, b_min, b_max)
        new(id_node, b_min, b_max, false, nothing)
    end
end

mutable struct Tree
    N::Int
    node::Vector{Node}
    function Tree(b_min, b_max)
        id_node = 1
        new(1, [Node(id_node, b_min, b_max)])
    end
end

function split!(tree::Tree, n_num::Int)
    # edit node
    leaves = [1, 2, 3, 4] .+ tree.N
    tree.node[n_num].hasLeaf = true
    tree.node[n_num].leaf = leaves

    # edit tree
    b_min = tree.node[n_num].b_min
    b_max = tree.node[n_num].b_max
    dif = b_max - b_min
    dx = [dif[1]*0.5, 0]
    dy = [0, dif[2]*0.5]
    b_center = b_min .+ dx .+ dy

    for i in 0:1, j in 0:1#, k in 0:1
        add = (i%2)*dx + (j%2)*dy 
        node_new = Node(tree.N+1, b_min+add, b_center+add)
        push!(tree.node, node_new)
        tree.N += 1
    end
end

function needSplitting(node::Node, f)
    v_lst = bound2vert(node.b_min, node.b_max)
    f_lst = [f(v) for v in v_lst]

    data = [f_lst[1] f_lst[4]; f_lst[2] f_lst[3]]
    itp = interpolate(data, BSpline(Linear()))

    center_itp = itp(1.5, 1.5) # note: interp starts from 1 
    center_real = f(0.5*(node.b_max + node.b_min))
    error = abs(center_itp - center_real)

    return ~(error<1.0) 
end

function auto_split!(tree::Tree, f) # recursive way
    # in the laef, we must add itp
    function recursion(n_node::Int)
        if needSplitting(tree.node[n_node], f)
            split!(tree, n_node)
            for id in tree.node[n_node].leaf
                recursion(id)
            end
        end
    end
    recursion(1)
    println("finish autosplit")
end

function show(tree::Tree)
    function recursion(n_node::Int)
        node = tree.node[n_node]
        if node.hasLeaf
            for id in tree.node[n_node].leaf
                recursion(id)
            end
        else
            v_lst = bound2vert(node.b_min, node.b_max)
            for n = 1:4
                if n!=4
                    x = [v_lst[n][1], v_lst[n+1][1]]
                    y = [v_lst[n][2], v_lst[n+1][2]]
                else
                    x = [v_lst[4][1], v_lst[1][1]]
                    y = [v_lst[4][2], v_lst[1][2]]
                end
                PyPlot.plot(x, y, "r-")
            end
        end
    end
    recursion(1)
    println("finish show")
end

f(p) = -norm(p)^2
t = Tree([-100, -100], [100, 100])
auto_split!(t, f)
show(t)
#needSplitting(t.node[1], f)
#split!(t, 1)

    

