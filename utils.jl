function bound2vert(b_min, b_max)
    dif = b_max - b_min
    dx = [dif[1], 0]
    dy = [0, dif[2]]

    v1 = b_min
    v2 = b_min + dx
    v3 = b_min + dx + dy
    v4 = b_min + dy
    v_lst = [v1, v2, v3, v4]
    return v_lst
end

function show(node::Node; color=:r)
    v_lst = bound2vert(node.b_min, node.b_max)
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

function itr()
    error("under construction")
    idx = 0
    function closure()
        add_x = (mod(idx, 2^1)==1)
        add_y = (mod(idx, 2^0)==1)
        idx += 1
    end
    return closure()
end


function whereami(q, b_min, b_max)
    ndim = length(q)
    idx = 1
    for n in 1:ndim
        tmp = (q[n] - b_min[n])/(b_max[n] - b_min[n])
        idx += (tmp>0.5)*2^(n-1)
    end
    return idx
end
