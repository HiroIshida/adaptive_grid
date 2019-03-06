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

function build_kdtree(vertex_lst)
    ndim = length(vertex_lst[1])
    N_vert = length(vertex_lst)
    
    # convert vector of vertex to 2dim array (mat) 
    vert_mat = zeros(ndim, N_vert)
    for n in 1:N_vert
        if ndim == 2
            vert_mat[:, n] = [vertex_lst[n][1], vertex_lst[n][2]]
        elseif ndim == 3
            vert_mat[:, n] = [vertex_lst[n][1], vertex_lst[n][2], vertex_lst[n][3]]
        end
    end
    kdtree = KDTree(vert_mat, leafsize=20)
    return kdtree
end

function vertex_reduction(vertex_lst)
    N_vert = length(vertex_lst)
    kdtree = build_kdtree(vertex_lst)
    id_lst = [i for i in 1:N_vert]

    map = [-1 for i in 1:N_vert] # -1 represetnts unvisited
    valid_ids = Int[] # thd ids extracted without depulication
    
    ε = 1e-4
    id_new = 1
    while(length(id_lst)>0)

        id = id_lst[1] # pop
        push!(valid_ids, id)

        # just display
        # just display
        percentage = id/N_vert*100
        #println(percentage)
        # just display
        # just display
        
        id_lst = setdiff(id_lst, id)
        map[id] = id_new
        id_depuli_lst = inrange(kdtree, vertex_lst[id], ε, true)
        for id_depuli in id_depuli_lst
            map[id_depuli] = id_new
            id_lst = setdiff(id_lst, id_depuli)
        end
        id_new += 1
    end

    return map, valid_ids
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

function grid_points(N_grid, b_min=[0, 0, 0], b_max=[1, 1, 1])
    # generate points at regular intervals in a cell 
    # return N_grid^3 points
    ndim = length(b_min)
    points = Vector{Float64}[]
    dx = (b_max - b_min)/N_grid
    if ndim == 2 # TODO
        for i in 0:N_grid-1, j in 0:N_grid-1
            p = ([i, j] .+ 0.5).*dx + b_min
            push!(points, p)
        end
    elseif ndim == 3
        for i in 0:N_grid-1, j in 0:N_grid-1, k in 0:N_grid-1
            p = ([i, j, k] .+ 0.5).*dx + b_min
            push!(points, p)
        end
    end
    return points
end

