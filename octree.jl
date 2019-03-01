using LinearAlgebra
using Interpolations
using SpecialFunctions
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
    id_child::Union{Vector{Int}, Nothing}
    itp # type?
    function Node(id_node, b_min, b_max)
        new(id_node, b_min, b_max, nothing, nothing)
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

function split!(tree::Tree, node::Node)
    # edit node
    leaves = [1, 2, 3, 4] .+ tree.N
    node.id_child = leaves

    # edit tree
    b_min = node.b_min
    b_max = node.b_max
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
                p_modif = (p - b_min)./(b_max - b_min) + 1
                return itp(p_modif[1], p_modif[2])
            end

        end
    end
    recursion(tree.node[1])
    println("finish autosplit")
end

function show(tree::Tree)
    function recursion(node::Node)
        if node.id_child!=nothing
            for id in node.id_child
                recursion(tree.node[id])
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
    recursion(tree.node[1])
    println("finish show")
end

function search_idx(tree::Tree, q)

    function recursion(n_node::Int)

    end

end



sigma = 8
f(x) = 0.5*(1 + erf((-norm(x)+40)/sqrt(2*sigma^2)))
t = Tree([-100, -100], [100, 100])
auto_split!(t, f)
show(t)

#=
function main1()
    sigma = 6
    f(x) = 1/sqrt(2*pi*sigma^2)*exp(-(norm(x)-30)^2/(2*sigma^2))
    x_start = 60
    v_lst = bound2vert([x_start, x_start], [x_start + 5, x_start + 5])
    f_lst = [f(v) for v in v_lst]
    data = [f_lst[1] f_lst[4]; f_lst[2] f_lst[3]]
    println(data)
    itp = interpolate(data, BSpline(Linear()))
    println(itp(1.5, 1.5))
    println(f([x_start + 2.5, x_start + 2.5]))
    v = itp(1.5, 1.5) - f([x_start + 2.5, x_start + 2.5])
    println(v)
end
main1()
=#
                  





    

