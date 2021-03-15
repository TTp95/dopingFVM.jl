"""
Boundary{T<:AbstractFloat, N <: Signed} <: TypeDopingFMV

"""
mutable struct BoundStructured{T<:AbstractFloat, N <: Signed} <: TypeDopingFMV
    kind::N
    cord::Char
    a::T
    b::T
    c::T
    d::T
    e::T
    eval::T
end
Integer
"""
Boundary{T<:AbstractFloat, N <: Signed} <: TypeDopingFMV

"""
struct BoundStructuredImmutable{T<:AbstractFloat, N <: Signed} <: TypeDopingFMV
    kind::N
    cord::Char
    a::T
    b::T
    c::T
    d::T
    e::T
    eval::T
end
