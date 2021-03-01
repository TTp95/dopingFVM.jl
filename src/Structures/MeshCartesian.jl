"""
CSMesh1D{T<:AbstractFloat, N <: Signed} <: CartesianStructuredMesh

"""
struct CSMesh1D{T<:AbstractFloat, N <: Signed} <: MeshCartesianStructured
    x::Array{T,1}
    dx::Array{T,1}
    vol::Array{T,1}
    l1::N
end

"""
CSMesh2D{T<:AbstractFloat, N <: Signed} <: CartesianStructuredMesh

"""
struct CSMesh2D{T<:AbstractFloat, N <: Signed} <: MeshCartesianStructured
    x::Array{T,1}
    y::Array{T,1}
    dx::Array{T,1}
    dy::Array{T,1}
    vol::Array{T,2}
    l1::N
    m1::N
end

"""
CSMesh3D{T<:AbstractFloat, N <: Signed} <: CartesianStructuredMesh

"""
struct CSMesh3D{T<:AbstractFloat, N <: Signed} <: MeshCartesianStructured
    x::Array{T,1}
    y::Array{T,1}
    z::Array{T,1}
    dx::Array{T,1}
    dy::Array{T,1}
    dz::Array{T,1}
    vol::Array{T,3}
    l1::N
    m1::N
    n1::N
end
