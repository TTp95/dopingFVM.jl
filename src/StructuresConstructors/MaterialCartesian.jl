"""

"""
function create_Material end

function create_Material(mesh::UnionCSMesh1D; T::Type{<:AbstractFloat}=Float64)
    return CSMaterial1D{T}(
        zeros(T, mesh.l1),
        zeros(T, mesh.l1),
        )
end

function create_Material(mesh::UnionCSMesh2D; T::Type{<:AbstractFloat}=Float64)
    return CSMaterial2D{T}(
        zeros(T, mesh.l1, mesh.m1),
        zeros(T, mesh.l1, mesh.m1),
        )
end

function create_Material(mesh::UnionCSMesh3D; T::Type{<:AbstractFloat}=Float64)
    return CSMaterial3D{T}(
        zeros(T, mesh.l1, mesh.m1, mesh.n1),
        zeros(T, mesh.l1, mesh.m1, mesh.n1),
        )
end
