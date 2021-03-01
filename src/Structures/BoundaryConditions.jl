"""
Boundary{T<:AbstractFloat, N <: Signed} <: BoundStructured

"""
struct BoundStructured{T<:AbstractFloat, N <: Signed} <: TypeDopingFMV
    kind::N
    cord::Char
    a::T
    b::T
    c::T
end
