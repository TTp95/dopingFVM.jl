"""
create_BoundStructured(;T=Float64, N=Int64)

"""
function create_BoundStructured(;
    T::Type{<:AbstractFloat}=Float64,
    N::Type{<:Signed}=Int64
)
    return BoundStructured{T, N}(0, 'c', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
end
