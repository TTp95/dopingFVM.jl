"""

"""
function create_CSMaterial end

function create_CSMaterial(mesh::UnionCSMesh1D; T::Type{<:AbstractFloat}=Float64)
    return CSMaterial1D{T}(
        zeros(T, mesh.l1),
        zeros(T, mesh.l1),
        )
end

function create_CSMaterial(mesh::UnionCSMesh2D; T::Type{<:AbstractFloat}=Float64)
    return CSMaterial2D{T}(
        zeros(T, mesh.l1, mesh.m1),
        zeros(T, mesh.l1, mesh.m1),
        )
end

function create_CSMaterial(mesh::UnionCSMesh3D; T::Type{<:AbstractFloat}=Float64)
    return CSMaterial3D{T}(
        zeros(T, mesh.l1, mesh.m1, mesh.n1),
        zeros(T, mesh.l1, mesh.m1, mesh.n1),
        )
end

function create_CSMaterial(; T::Type{<:AbstractFloat}=Float64, mutable::Bool=true)
    if mutable
        return CSMaterialConstant{T}(0.0, 0.0)
    elseif !mutable
        return CSMaterialConstantImmutable{T}(0.0, 0.0)
    end
end

function create_CSMaterial(
    ρ::AbstractFloat,
    Γ::AbstractFloat;
    T::Type{<:AbstractFloat}=Float64,
    mutable::Bool=true,
)
    if mutable
        return CSMaterialConstant{T}(ρ, Γ)
    elseif !mutable
        return CSMaterialConstantImmutable{T}(ρ, Γ)
    end
end
