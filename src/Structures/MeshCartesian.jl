"""
CSMesh1D{T<:AbstractFloat, N <: Signed} <: StructuredMesh

"""
mutable struct CSMesh1D{T<:AbstractFloat, N <: Signed} <: MeshStructured
    x::Array{T,1}
    dx::Array{T,1}
    vol::Array{T,1}
    l1::N
end

"""
CSMesh2D{T<:AbstractFloat, N <: Signed} <: MeshStructured

"""
mutable struct CSMesh2D{T<:AbstractFloat, N <: Signed} <: MeshStructured
    x::Array{T,1}
    y::Array{T,1}
    dx::Array{T,1}
    dy::Array{T,1}
    vol::Array{T,2}
    l1::N
    m1::N
end

"""
CSMesh3D{T<:AbstractFloat, N <: Signed} <: MeshStructured

"""
mutable struct CSMesh3D{T<:AbstractFloat, N <: Signed} <: MeshStructured
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

"""
CSMesh1D{T<:AbstractFloat, N <: Signed} <: MeshStructured

"""
struct CSMesh1DImmutable{T<:AbstractFloat, N <: Signed} <: MeshStructured
    x::Array{T,1}
    dx::Array{T,1}
    vol::Array{T,1}
    l1::N
end

"""
CSMesh2D{T<:AbstractFloat, N <: Signed} <: MeshStructured

"""
struct CSMesh2DImmutable{T<:AbstractFloat, N <: Signed} <: MeshStructured
    x::Array{T,1}
    y::Array{T,1}
    dx::Array{T,1}
    dy::Array{T,1}
    vol::Array{T,2}
    l1::N
    m1::N
end

"""
CSMesh3DImmutable{T<:AbstractFloat, N <: Signed} <: MeshStructured

"""
struct CSMesh3DImmutable{T<:AbstractFloat, N <: Signed} <: MeshStructured
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
