"""

"""
function discretize_convection end

@inline function discretize_convection(
    vel::CSVelocity1D,
    phi::CSPhi1D,
    bounds::Dict{String,BoundsStructured},
    material::Union{CSMaterial1D,UnionCSConstantMaterial},
    mesh::UnionCSMesh1D;
    Adiff::Union{
        SparseVector{<:AbstractFloat,<:Signed},
        SparseMatrixCSC{<:AbstractFloat,<:Signed},
        Array{<:AbstractFloat,2},
    } = [0.0 0.0; 0.0 0.0;],
    velocityU::Array{<:AbstractFloat,1} = vel.fValues.uFace,
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
    mthreads::Bool = false,
    sparrays::Bool = true,
    scheme::Signed = 5,
    interpolation::Signed = 1,
)
    if (scheme == 1)
        A, b = _discretize_convection_centralDifference_(
            vel,
            phi,
            bounds,
            material,
            mesh;
            velocityU = velocityU,
            T = T,
            N = N,
            mthreads = mthreads,
            sparrays = sparrays,
            interpolation = interpolation,
        )

    elseif (scheme == 2)
        A, b = _discretize_convection_upwind_(
            vel,
            phi,
            bounds,
            material,
            mesh;
            velocityU = velocityU,
            T = T,
            N = N,
            mthreads = mthreads,
            sparrays = sparrays,
            interpolation = interpolation,
        )

    elseif (scheme == 3)
        A, b = _discretize_convection_downwind_(
            vel,
            phi,
            bounds,
            material,
            mesh;
            velocityU = velocityU,
            T = T,
            N = N,
            mthreads = mthreads,
            sparrays = sparrays,
            interpolation = interpolation,
        )

    elseif (scheme == 4)
        A, b = _discretize_convection_hybrid_(
            vel,
            phi,
            bounds,
            material,
            mesh;
            Adiff = Adiff,
            velocityU = velocityU,
            T = T,
            N = N,
            mthreads = mthreads,
            sparrays = sparrays,
            interpolation = interpolation,
        )

    elseif (scheme == 5)
        A, b = _discretize_convection_powerlaw_(
            vel,
            phi,
            bounds,
            material,
            mesh;
            Adiff = Adiff,
            velocityU = velocityU,
            T = T,
            N = N,
            mthreads = mthreads,
            sparrays = sparrays,
            interpolation = interpolation,
        )

    elseif (scheme == 6)
        A, b = _discretize_convection_secondorderupwind_(
            vel,
            phi,
            bounds,
            material,
            mesh;
            velocityU = velocityU,
            T = T,
            N = N,
            mthreads = mthreads,
            sparrays = sparrays,
            interpolation = interpolation,
        )

    else
        error("Diffusion scheme number $(scheme) unimplemented.")
    end

    return A, b
end

@inline function discretize_convection(
    vel::CSVelocity2D,
    phi::CSPhi2D,
    bounds::Dict{String,BoundsStructured},
    material::Union{CSMaterial2D,UnionCSConstantMaterial},
    mesh::UnionCSMesh2D;
    Adiff::Union{
        SparseVector{<:AbstractFloat,<:Signed},
        SparseMatrixCSC{<:AbstractFloat,<:Signed},
        Array{<:AbstractFloat,2},
    } = [0.0 0.0; 0.0 0.0;],
    velocityU::Array{<:AbstractFloat,2} = vel.fValues.uFace,
    velocityV::Array{<:AbstractFloat,2} = vel.fValues.vFace,
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
    mthreads::Bool = false,
    sparrays::Bool = true,
    scheme::Signed = 5,
    interpolation::Signed = 1,
)
    if (scheme == 1)
        A, b = _discretize_convection_centralDifference_(
            vel,
            phi,
            bounds,
            material,
            mesh;
            velocityU = velocityU,
            velocityV = velocityV,
            T = T,
            N = N,
            mthreads = mthreads,
            sparrays = sparrays,
            interpolation = interpolation,
        )

    elseif (scheme == 2)
        A, b = _discretize_convection_upwind_(
            vel,
            phi,
            bounds,
            material,
            mesh;
            velocityU = velocityU,
            velocityV = velocityV,
            T = T,
            N = N,
            mthreads = mthreads,
            sparrays = sparrays,
            interpolation = interpolation,
        )

    elseif (scheme == 3)
        A, b = _discretize_convection_downwind_(
            vel,
            phi,
            bounds,
            material,
            mesh;
            velocityU = velocityU,
            velocityV = velocityV,
            T = T,
            N = N,
            mthreads = mthreads,
            sparrays = sparrays,
            interpolation = interpolation,
        )

    elseif (scheme == 4)
        A, b = _discretize_convection_hybrid_(
            vel,
            phi,
            bounds,
            material,
            mesh;
            Adiff = Adiff,
            velocityU = velocityU,
            velocityV = velocityV,
            T = T,
            N = N,
            mthreads = mthreads,
            sparrays = sparrays,
            interpolation = interpolation,
        )

    elseif (scheme == 5)
        A, b = _discretize_convection_powerlaw_(
            vel,
            phi,
            bounds,
            material,
            mesh;
            Adiff = Adiff,
            velocityU = velocityU,
            velocityV = velocityV,
            T = T,
            N = N,
            mthreads = mthreads,
            sparrays = sparrays,
            interpolation = interpolation,
        )

    elseif (scheme == 6)
        A, b = _discretize_convection_secondorderupwind_(
            vel,
            phi,
            bounds,
            material,
            mesh;
            velocityU = velocityU,
            velocityV = velocityV,
            T = T,
            N = N,
            mthreads = mthreads,
            sparrays = sparrays,
            interpolation = interpolation,
        )

    else
        error("Diffusion scheme number $(scheme) unimplemented.")
    end

    return A, b
end

@inline function discretize_convection(
    vel::CSVelocity3D,
    phi::CSPhi3D,
    bounds::Dict{String,BoundsStructured},
    material::Union{CSMaterial3D,UnionCSConstantMaterial},
    mesh::UnionCSMesh3D;
    Adiff::Union{
        SparseVector{<:AbstractFloat,<:Signed},
        SparseMatrixCSC{<:AbstractFloat,<:Signed},
        Array{<:AbstractFloat,2},
    } = [0.0 0.0; 0.0 0.0;],
    velocityU::Array{<:AbstractFloat,3} = vel.fValues.uFace,
    velocityV::Array{<:AbstractFloat,3} = vel.fValues.vFace,
    velocityW::Array{<:AbstractFloat,3} = vel.fValues.wFace,
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
    mthreads::Bool = false,
    sparrays::Bool = true,
    scheme::Signed = 5,
    interpolation::Signed = 1,
)
    if (scheme == 1)
        A, b = _discretize_convection_centralDifference_(
            vel,
            phi,
            bounds,
            material,
            mesh;
            velocityU = velocityU,
            velocityV = velocityV,
            velocityW = velocityW,
            T = T,
            N = N,
            mthreads = mthreads,
            sparrays = sparrays,
            interpolation = interpolation,
        )

    elseif (scheme == 2)
        A, b = _discretize_convection_upwind_(
            vel,
            phi,
            bounds,
            material,
            mesh;
            velocityU = velocityU,
            velocityV = velocityV,
            velocityW = velocityW,
            T = T,
            N = N,
            mthreads = mthreads,
            sparrays = sparrays,
            interpolation = interpolation,
        )

    elseif (scheme == 3)
        A, b = _discretize_convection_downwind_(
            vel,
            phi,
            bounds,
            material,
            mesh;
            velocityU = velocityU,
            velocityV = velocityV,
            velocityW = velocityW,
            T = T,
            N = N,
            mthreads = mthreads,
            sparrays = sparrays,
            interpolation = interpolation,
        )

    elseif (scheme == 4)
        A, b = _discretize_convection_hybrid_(
            vel,
            phi,
            bounds,
            material,
            mesh;
            Adiff = Adiff,
            velocityU = velocityU,
            velocityV = velocityV,
            velocityW = velocityW,
            T = T,
            N = N,
            mthreads = mthreads,
            sparrays = sparrays,
            interpolation = interpolation,
        )

    elseif (scheme == 5)
        A, b = _discretize_convection_powerlaw_(
            vel,
            phi,
            bounds,
            material,
            mesh;
            Adiff = Adiff,
            velocityU = velocityU,
            velocityV = velocityV,
            velocityW = velocityW,
            T = T,
            N = N,
            mthreads = mthreads,
            sparrays = sparrays,
            interpolation = interpolation,
        )

    elseif (scheme == 6)
        A, b = _discretize_convection_secondorderupwind_(
            vel,
            phi,
            bounds,
            material,
            mesh;
            velocityU = velocityU,
            velocityV = velocityV,
            velocityW = velocityW,
            T = T,
            N = N,
            mthreads = mthreads,
            sparrays = sparrays,
            interpolation = interpolation,
        )

    else
        error("Diffusion scheme number $(scheme) unimplemented.")
    end

    return A, b
end
