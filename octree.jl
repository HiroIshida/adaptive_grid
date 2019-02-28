using LinearAlgebra
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
    function Tree(f, b_min, b_max)
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

f(p) = norm(p)
t = Tree(f, [-10, -10], [10, 10])
split!(t, 1)
    

