"""

"""
function discretize_SIMPLE_PressureCorrection end

function discretize_SIMPLE_PressureCorrection(
    AU::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,2},
    },
    velocityU::Array{<:AbstractFloat,1},
    phi::CSPhi1D,
    material::Union{CSMaterial1D,UnionCSConstantMaterial},
    mesh::UnionCSMesh1D;
    T::Type{<:AbstractFloat} = Float64,
    threads::Bool = false,
    sparse::Bool = true,
    interpolation::Signed = 1,
)
    n_equations = maximum_globalIndex(phi)

    if sparse
        A = spzeros(T, n_equations, n_equations)
    elseif !sparse
        A = zeros(T, n_equations, n_equations)
    end

    b = zeros(T, n_equations)

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            if phi.onoff[i]
                id = phi.gIndex[i]
                @inbounds A[id,:], b[id] = _SIMPLE_PressureCorrection_Coefficients_(
                    AU,
                    i,
                    velocityU,
                    phi,
                    material,
                    mesh;
                    T = T,
                    sparse = sparse,
                    interpolation = interpolation,
                )
            end
        end

    elseif !threads
        for i in 1:mesh.l1
            if phi.onoff[i]
                id = phi.gIndex[i]
                @inbounds A[id,:], b[id] = _SIMPLE_PressureCorrection_Coefficients_(
                    AU,
                    i,
                    velocityU,
                    phi,
                    material,
                    mesh;
                    T = T,
                    sparse = sparse,
                    interpolation = interpolation,
                )
            end
        end

    end

    if sparse
        A = dropzeros(A)
    end

    return A, b
end

function discretize_SIMPLE_PressureCorrection(
    AU::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,2},
    },
    AV::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,2},
    },
    velocityU::Array{<:AbstractFloat,2},
    velocityV::Array{<:AbstractFloat,2},
    phi::CSPhi2D,
    material::Union{CSMaterial2D,UnionCSConstantMaterial},
    mesh::UnionCSMesh2D;
    T::Type{<:AbstractFloat} = Float64,
    threads::Bool = false,
    sparse::Bool = true,
    scheme::Signed = 2,
    interpolation::Signed = 1,
)
    n_equations = maximum_globalIndex(phi)

    if sparse
        A = spzeros(T, n_equations, n_equations)
    elseif !sparse
        A = zeros(T, n_equations, n_equations)
    end

    b = zeros(T, n_equations)

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            for j in 1:mesh.m1
                if phi.onoff[i,j]
                    id = phi.gIndex[i,j]
                    @inbounds A[id,:], b[id] = _SIMPLE_PressureCorrection_Coefficients_(
                        AU,
                        AV,
                        i,
                        j,
                        velocityU,
                        velocityV,
                        phi,
                        material,
                        mesh;
                        T = T,
                        sparse = sparse,
                        interpolation = interpolation,
                    )
                end
            end
        end

    elseif !threads
        for i in 1:mesh.l1
            for j in 1:mesh.m1
                if phi.onoff[i,j]
                    id = phi.gIndex[i,j]
                    @inbounds A[id,:], b[id] = _SIMPLE_PressureCorrection_Coefficients_(
                        AU,
                        AV,
                        i,
                        j,
                        velocityU,
                        velocityV,
                        phi,
                        material,
                        mesh;
                        T = T,
                        sparse = sparse,
                        interpolation = interpolation,
                    )
                end
            end
        end

    end

    if sparse
        A = dropzeros(A)
    end

    return A, b
end

function discretize_SIMPLE_PressureCorrection(
    AU::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,2},
    },
    AV::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,2},
    },
    AW::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,2},
    },
    velocityU::Array{<:AbstractFloat,3},
    velocityV::Array{<:AbstractFloat,3},
    velocityW::Array{<:AbstractFloat,3},
    phi::CSPhi3D,
    material::Union{CSMaterial3D,UnionCSConstantMaterial},
    mesh::UnionCSMesh3D;
    T::Type{<:AbstractFloat} = Float64,
    threads::Bool = false,
    sparse::Bool = true,
    scheme::Signed = 2,
    interpolation::Signed = 1,
)
    n_equations = maximum_globalIndex(phi)

    if sparse
        A = spzeros(T, n_equations, n_equations)
    elseif !sparse
        A = zeros(T, n_equations, n_equations)
    end

    b = zeros(T, n_equations)

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            for j in 1:mesh.m1
                for k in 1:mesh.n1
                    if phi.onoff[i,j,k]
                        id = phi.gIndex[i,j,k]
                        @inbounds A[id,:], b[id] = _SIMPLE_PressureCorrection_Coefficients_(
                            AU,
                            AV,
                            AW,
                            i,
                            j,
                            k,
                            velocityU,
                            velocityV,
                            velocityW,
                            phi,
                            material,
                            mesh;
                            T = T,
                            sparse = sparse,
                            interpolation = interpolation,
                        )
                    end
                end
            end
        end

    elseif !threads
        for i in 1:mesh.l1
            for j in 1:mesh.m1
                for k in 1:mesh.n1
                    if phi.onoff[i,j,k]
                        id = phi.gIndex[i,j,k]
                        @inbounds A[id,:], b[id] = _SIMPLE_PressureCorrection_Coefficients_(
                            AU,
                            AV,
                            AW,
                            i,
                            j,
                            k,
                            velocityU,
                            velocityV,
                            velocityW,
                            phi,
                            material,
                            mesh;
                            T = T,
                            sparse = sparse,
                            interpolation = interpolation,
                        )
                    end
                end
            end
        end

    end

    if sparse
        A = dropzeros(A)
    end

    return A, b
end
