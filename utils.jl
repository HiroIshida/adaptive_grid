function bound2vert(b_min, b_max)
    ndim = length(b_min)
    dx = bound2dx(b_min, b_max)
    v_lst = Vector{Float64}[]

    for i in 0:(2^ndim-1)
        add = [0.0 for n in 1:ndim]
        for dim in 1:ndim
            if mod(div(i, 2^(dim-1)), 2) == 1
                add += dx[dim]
            end
        end
        push!(v_lst, b_min + add)

    end

    #=
    v1 = b_min
    v2 = b_min + dx
    v3 = b_min + dx + dy
    v4 = b_min + dy
    v_lst = [v1, v2, v3, v4]
    =#
    return v_lst
end

function bound2dx(b_min, b_max)
    ndim = length(b_min)
    dif = b_max - b_min
    dx = Vector{Float64}[]
    for i in 1:ndim
        dx_ = [0.0 for n=1:ndim]
        dx_[i] = dif[i]
        push!(dx, dx_)
    end
    return dx
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
