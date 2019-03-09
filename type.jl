using StaticArrays
struct A{N, T}
    v::SVector{N, T}
    function A(v)
        new{length(v), Float64}(v)
    end
end

b = A(SVector(4 ,2, 2))

