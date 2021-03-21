"""
create_BoundsStructured(;T=Float64, N=Int64)

"""
function create_BoundsStructured end

function create_BoundsStructured(
    ;
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
)
    return BoundsStructured{T, N}(0, 'c', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
end

function create_BoundsStructured(
    kind::Signed,
    cord::Char,
    ρ::AbstractFloat,
    Γ::AbstractFloat,
    a::AbstractFloat,
    b::AbstractFloat,
    c::AbstractFloat,
    eval::AbstractFloat;
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
)
    return BoundsStructured{T, N}(kind, cord, ρ, Γ, a, b, c, eval)
end
