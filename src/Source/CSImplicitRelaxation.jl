"""

"""
function implicit_relaxation! end

function implicit_relaxation!(
    A::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,2},
    },
    b::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,1},
    },
    phi::CSPhi1D,
    mesh::UnionCSMesh1D;
    relax = 1.0,
    T::Type{<:AbstractFloat} = Float64,
    threads::Bool = false,
)
    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            if phi.onoff[i]
                id = phi.gIndex[i]
                @inbounds b[id]  += ((1.0 - relax) / relax) * A[id,id] * phi.iter[i]
                @inbounds A[id,id]  = A[id,id] / relax
            end
        end
    elseif !threads
        Base.Threads.@threads for i in 1:mesh.l1
            if phi.onoff[i]
                id = phi.gIndex[i]
                @inbounds b[id]  += ((1.0 - relax) / relax) * A[id,id] * phi.iter[i]
                @inbounds A[id,id]  = A[id,id] / relax
            end
        end
    end

    return nothing
end

function implicit_relaxation!(
    A::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,2},
    },
    b::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,1},
    },
    phi::CSPhi2D,
    mesh::UnionCSMesh2D;
    relax = 1.0,
    T::Type{<:AbstractFloat} = Float64,
    threads::Bool = false,
)
    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            for j in 1:mesh.m1
                if phi.onoff[i,j]
                    id = phi.gIndex[i,j]
                    @inbounds b[id]  += ((1.0 - relax) / relax) * A[id,id] * phi.iter[i,j]
                    @inbounds A[id,id]  = A[id,id] / relax
                end
            end
        end
    elseif !threads
        for i in 1:mesh.l1
            for j in 1:mesh.m1
                if phi.onoff[i,j]
                    id = phi.gIndex[i,j]
                    @inbounds b[id]  += ((1.0 - relax) / relax) * A[id,id] * phi.iter[i,j]
                    @inbounds A[id,id]  = A[id,id] / relax
                end
            end
        end
    end

    return nothing
end

function implicit_relaxation!(
    A::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,2},
    },
    b::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,1},
    },
    phi::CSPhi3D,
    mesh::UnionCSMesh3D;
    relax = 1.0,
    T::Type{<:AbstractFloat} = Float64,
    threads::Bool = false,
)
    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            for j in 1:mesh.m1
                for k in 1:mesh.n1
                    if phi.onoff[i,j,k]
                        id = phi.gIndex[i,j,k]
                        @inbounds b[id]  += ((1.0 - relax) / relax) * A[id,id] * phi.iter[i,j,k]
                        @inbounds A[id,id]  = A[id,id] / relax
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
                        @inbounds b[id]  += ((1.0 - relax) / relax) * A[id,id] * phi.iter[i,j,k]
                        @inbounds A[id,id]  = A[id,id] / relax
                    end
                end
            end
        end
    end

    return nothing
end

"""

"""
function implicit_relaxation end

function implicit_relaxation(
    A::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,2},
    },
    phi::CSPhi1D,
    mesh::UnionCSMesh1D;
    relax = 1.0,
    T::Type{<:AbstractFloat} = Float64,
    threads::Bool = false,
    sparse::Bool = true,
)
    n_equations = maximum_globalIndex(phi)

    if sparse
        Ar = spzeros(T, n_equations, n_equations)
    elseif !sparse
        Ar = zeros(T, n_equations, n_equations)
    end

    br = zeros(T, n_equations)

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            if phi.onoff[i]
                id = phi.gIndex[i]
                @inbounds br[id] = ((1.0 - relax) / relax) * A[id,id] * phi.iter[i]
                @inbounds Ar[id,id]  = (A[id,id] / relax) - A[id,id]
            end
        end
    elseif !threads
        Base.Threads.@threads for i in 1:mesh.l1
            if phi.onoff[i]
                id = phi.gIndex[i]
                @inbounds br[id] = ((1.0 - relax) / relax) * A[id,id] * phi.iter[i]
                @inbounds Ar[id,id]  = (A[id,id] / relax) - A[id,id]
            end
        end
    end

    return Ar, br
end

function implicit_relaxation(
    A::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,2},
    },
    b::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,1},
    },
    phi::CSPhi2D,
    mesh::UnionCSMesh2D;
    relax = 1.0,
    T::Type{<:AbstractFloat} = Float64,
    threads::Bool = false,
)
    n_equations = maximum_globalIndex(phi)

    if sparse
        Ar = spzeros(T, n_equations, n_equations)
    elseif !sparse
        Ar = zeros(T, n_equations, n_equations)
    end

    br = zeros(T, n_equations)

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            for j in 1:mesh.m1
                if phi.onoff[i,j]
                    id = phi.gIndex[i,j]
                    @inbounds br[id] = ((1.0 - relax) / relax) * A[id,id] * phi.iter[i,j]
                    @inbounds Ar[id,id]  = (A[id,id] / relax) - A[id,id]
                end
            end
        end
    elseif !threads
        for i in 1:mesh.l1
            for j in 1:mesh.m1
                if phi.onoff[i,j]
                    id = phi.gIndex[i,j]
                    @inbounds br[id] = ((1.0 - relax) / relax) * A[id,id] * phi.iter[i,j]
                    @inbounds Ar[id,id]  = (A[id,id] / relax) - A[id,id]
                end
            end
        end
    end

    return Ar, br
end

function implicit_relaxation(
    A::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,2},
    },
    b::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,1},
    },
    phi::CSPhi3D,
    mesh::UnionCSMesh3D;
    relax = 1.0,
    T::Type{<:AbstractFloat} = Float64,
    threads::Bool = false,
)
    n_equations = maximum_globalIndex(phi)

    if sparse
        Ar = spzeros(T, n_equations, n_equations)
    elseif !sparse
        Ar = zeros(T, n_equations, n_equations)
    end

    br = zeros(T, n_equations)

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            for j in 1:mesh.m1
                for k in 1:mesh.n1
                    if phi.onoff[i,j,k]
                        id = phi.gIndex[i,j,k]
                        @inbounds br[id] = ((1.0 - relax) / relax) * A[id,id] * phi.iter[i,j,k]
                        @inbounds Ar[id,id] = (A[id,id] / relax) - A[id,id]
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
                        @inbounds br[id] = ((1.0 - relax) / relax) * A[id,id] * phi.iter[i,j,k]
                        @inbounds Ar[id,id] = (A[id,id] / relax) - A[id,id]
                    end
                end
            end
        end
    end

    return Ar, br
end
