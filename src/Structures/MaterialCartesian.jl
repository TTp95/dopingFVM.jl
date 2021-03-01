"""
CSMaterial1D{T<:AbstractFloat, N <: Signed} <: MaterialCartesianStructured

"""
struct CSMaterial1D{T<:AbstractFloat, N <: Signed} <: MaterialCartesianStructured
    ρ::Array{T,1}
    μ::Array{T,1}
    κ::Array{T,1}
    cp::Array{T,1}
end

"""
CSMaterialConstant1D{T<:AbstractFloat, N <: Signed} <: MaterialCartesianStructured

"""
struct CSMaterialConstant1D{T<:AbstractFloat, N <: Signed} <: MaterialCartesianStructured
    ρ::T
    μ::T
    κ::T
    cp::T
end

"""
CSMaterial2D{T<:AbstractFloat, N <: Signed} <: MaterialCartesianStructured

"""
struct CSMaterial2D{T<:AbstractFloat, N <: Signed} <: MaterialCartesianStructured
    ρ::Array{T,2}
    μ::Array{T,2}
    κ::Array{T,2}
    cp::Array{T,2}
end

"""
CSMaterialConstant2D{T<:AbstractFloat, N <: Signed} <: MaterialCartesianStructured

"""
struct CSMaterialConstant2D{T<:AbstractFloat, N <: Signed} <: MaterialCartesianStructured
    ρ::T
    μ::T
    κ::T
    cp::T
end

"""
CSMaterial3D{T<:AbstractFloat, N <: Signed} <: MaterialCartesianStructured

"""
struct CSMaterial3D{T<:AbstractFloat, N <: Signed} <: MaterialCartesianStructured
    ρ::Array{T,3}
    μ::Array{T,3}
    κ::Array{T,3}
    cp::Array{T,3}
end

"""
CSMaterialConstant3D{T<:AbstractFloat, N <: Signed} <: MaterialCartesianStructured

"""
struct CSMaterialConstant3D{T<:AbstractFloat, N <: Signed} <: MaterialCartesianStructured
    ρ::T
    μ::T
    κ::T
    cp::T
end
