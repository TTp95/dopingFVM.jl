"""
CSMaterial1D{T<:AbstractFloat} <: MaterialCartesianStructured

"""
struct CSMaterial1D{T<:AbstractFloat} <: MaterialCartesianStructured
    ρ::Array{T,1}
    Γ::Array{T,1}
end

"""
CSMaterial2D{T<:AbstractFloat} <: MaterialCartesianStructured

"""
struct CSMaterial2D{T<:AbstractFloat} <: MaterialCartesianStructured
    ρ::Array{T,2}
    Γ::Array{T,2}
end

"""
CSMaterial3D{T<:AbstractFloat} <: MaterialCartesianStructured

"""
struct CSMaterial3D{T<:AbstractFloat} <: MaterialCartesianStructured
    ρ::Array{T,3}
    Γ::Array{T,3}
end

"""
CSMaterialConstant{T<:AbstractFloat} <: MaterialCartesianStructured

"""
mutable struct CSMaterialConstant{T<:AbstractFloat} <: MaterialCartesianStructured
    ρ::T
    Γ::T
end

"""
CSMaterialConstant{T<:AbstractFloat} <: MaterialCartesianStructured

"""
struct CSMaterialConstantImmutable{T<:AbstractFloat} <: MaterialCartesianStructured
    ρ::T
    Γ::T
end
