include("octree.jl")
using PyPlot
import Base: show

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
