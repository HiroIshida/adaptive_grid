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
    return v_lst
end

function itr(ndim)
    i = 0
    return function closure()
        indicater = Int[]
        for dim in 1:ndim
            push!(indicater, mod(div(i, 2^(dim-1)), 2))
        end
        i += 1
        return (i<2^ndim+1 ? indicater : nothing)
    end

end

function form_data_cubic(f_lst, ndim)
    # TODO fix this by using @eval
    # just for now
    if ndim == 2
        data = zeros(2, 2)
    elseif ndim == 3
        data = zeros(2, 2, 2)
    end

    # TODO dirty dirty dirty dirty
    it = itr(ndim)
    for i = 1:2^ndim
        idc = it()
        if ndim == 2
            data[idc[1]+1, idc[2]+1] = f_lst[i]
        elseif ndim == 3
            data[idc[1]+1, idc[2]+1, idc[3]+1] = f_lst[i]
        end
    end
    return data
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

function grid_points(N_grid, b_min, b_max)
    # generate points at regular intervals in a cell 
    # return N_grid^3 points
    points = Vector{Float64}[]
    dx = (b_max - b_min)/N_grid
    for i in 0:N_grid-1, j in 0:N_grid-1, k in 0:N_grid-1
        p = ([i, j, k] .+ 0.5).*dx + b_min
        push!(points, p)
    end
    return points
end
