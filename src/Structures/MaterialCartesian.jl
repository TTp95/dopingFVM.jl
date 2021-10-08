"""
CSMaterial1D{T<:AbstractFloat} <: MaterialStructured

"""
struct CSMaterial1D{T<:AbstractFloat} <: MaterialStructured
    ρ::Array{T,1}
    Γ::Array{T,1}
end

"""
CSMaterial2D{T<:AbstractFloat} <: MaterialStructured

"""
struct CSMaterial2D{T<:AbstractFloat} <: MaterialStructured
    ρ::Array{T,2}
    Γ::Array{T,2}
end

"""
CSMaterial3D{T<:AbstractFloat} <: MaterialStructured

"""
struct CSMaterial3D{T<:AbstractFloat} <: MaterialStructured
    ρ::Array{T,3}
    Γ::Array{T,3}
end

