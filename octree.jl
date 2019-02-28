using LinearAlgebra
using Interpolations
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
    dif = node.b_max - node.b_min
    dx = [dif[1]*0.5, 0]
    dy = [0, dif[2]*0.5]
    
    v1 = node.b_min; f1 = f(v1)
    v2 = node.b_min + dx; f2 = f(v2)
    v3 = node.b_min + dx + dy; f3 = f(v3)
    v4 = node.b_min + dy; f4 = f(v4)

    data = [f1 f4; f2 f3]
    itp = interpolate(data, BSpline(Linear()))

    center_itp = itp(1.5, 1.5) # note: interp starts from 1 
    center_real = f(0.5*(node.b_max - node.b_min))

    return ~(abs(center_itp - center_real)<0.1)
end

function split_grid(tree::Tree) # recursive way
    function recursion(node::Node)
        if needSplitting
        end
    end
end

f(p) = -norm(p)^2
t = Tree([-2, -2], [2, 2])
needSplitting(t.node[1], f)
#split!(t, 1)

    

