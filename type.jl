using StaticArrays
struct A{N, T}
    v::SVector{N, T}
    function A(v)
        new{length(v), Float64}(v)
    end
end

function fuck(a::SVector{N, Float64}) where N
    println(N)
    return 1
end

N=2
fuck(SVector(2, 2.0)::SVector{N, Float64})

