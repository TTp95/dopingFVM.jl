"""
CSMaterial1D{T<:AbstractFloat} <: MaterialCartesianStructured

"""
struct CSMaterial1D{T<:AbstractFloat} <: MaterialCartesianStructured
    ρ::Array{T,1}
    μ::Array{T,1}
    κ::Array{T,1}
    cp::Array{T,1}
end

"""
CSMaterial2D{T<:AbstractFloat} <: MaterialCartesianStructured

"""
struct CSMaterial2D{T<:AbstractFloat} <: MaterialCartesianStructured
    ρ::Array{T,2}
    μ::Array{T,2}
    κ::Array{T,2}
    cp::Array{T,2}
end

"""
CSMaterial3D{T<:AbstractFloat} <: MaterialCartesianStructured

"""
struct CSMaterial3D{T<:AbstractFloat} <: MaterialCartesianStructured
    ρ::Array{T,3}
    μ::Array{T,3}
    κ::Array{T,3}
    cp::Array{T,3}
end

"""
CSMaterialConstant{T<:AbstractFloat} <: MaterialCartesianStructured

"""
struct CSMaterialConstant{T<:AbstractFloat} <: MaterialCartesianStructured
    ρ::T
    μ::T
    κ::T
    cp::T
end
