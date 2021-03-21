"""
Boundary{T<:AbstractFloat, N <: Signed} <: TypeDopingFMV

"""
mutable struct BoundsStructured{T <: AbstractFloat, N <: Signed} <: TypeDopingFMV
    kind::N
    cord::Char
    ρ::T
    Γ::T
    a::T
    b::T
    c::T
    eval::T
end
