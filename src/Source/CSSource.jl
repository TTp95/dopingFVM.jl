"""

"""
function discretize_source end

function discretize_source(
    phi::CSPhi1D,
    mesh::UnionCSMesh1D;
    T::Type{<:AbstractFloat} = Float64,
    threads::Bool = false,
    sparse::Bool = true,
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
                @inbounds A[id,id]  = -1.0 * phi.sourceP[i] * mesh.vol[i]
                @inbounds b[id]  = phi.sourceC[i] * mesh.vol[i]
            end
        end
    elseif !threads
        for i in 1:mesh.l1
            if phi.onoff[i]
                id = phi.gIndex[i]
                @inbounds A[id,id]  = -1.0 * phi.sourceP[i] * mesh.vol[i]
                @inbounds b[id]  = phi.sourceC[i] * mesh.vol[i]
            end
        end
    end

    return A, b
end

function discretize_source(
    phi::CSPhi2D,
    mesh::UnionCSMesh2D;
    T::Type{<:AbstractFloat} = Float64,
    threads::Bool = false,
    sparse::Bool = true,
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
                    @inbounds A[id,id]  = -1.0 * phi.sourceP[i,j] * mesh.vol[i,j]
                    @inbounds b[id]  = phi.sourceC[i,j] * mesh.vol[i,j]
                end
            end
        end
    elseif !threads
        for i in 1:mesh.l1
            for j in 1:mesh.m1
                if phi.onoff[i,j]
                    id = phi.gIndex[i,j]
                    @inbounds A[id,id]  = -1.0 * phi.sourceP[i,j] * mesh.vol[i,j]
                    @inbounds b[id]  = phi.sourceC[i,j] * mesh.vol[i,j]
                end
            end
        end
    end

    return A, b
end

function discretize_source(
    phi::CSPhi3D,
    mesh::UnionCSMesh3D;
    T::Type{<:AbstractFloat} = Float64,
    threads::Bool = false,
    sparse::Bool = true,
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
                        @inbounds A[id,id]  = -1.0 * phi.sourceP[i,j,k] * mesh.vol[i,j,k]
                        @inbounds b[id]  = phi.sourceC[i,j,k] * mesh.vol[i,j,k]
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
                        @inbounds A[id,id]  = -1.0 * phi.sourceP[i,j,k] * mesh.vol[i,j,k]
                        @inbounds b[id]  = phi.sourceC[i,j,k] * mesh.vol[i,j,k]
                    end
                end
            end
        end
    end

    return A, b
end
