"""

"""
function vector_to_phi! end

function vector_to_phi!(
    phi::CSPhi1D,
    mesh::UnionCSMesh1D;
    phisolution::Array{<:AbstractFloat,1} = phi.eval,
    T::Type{<:AbstractFloat} = Float64,
    threads = false,
)
    n_equations = maximum_globalIndex(phi)

    x = zeros(T, n_equations)

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            if phi.onoff[i]
                id = phi.gIndex[i]
                x[id] = phisolution[i]
            end
        end
    elseif !threads
        for i in 1:mesh.l1
            if phi.onoff[i]
                id = phi.gIndex[i]
                x[id] = phisolution[i]
            end
        end
    end

    return x
end

function vector_to_phi!(
    phi::CSPhi2D,
    mesh::UnionCSMesh2D;
    phisolution::Array{<:AbstractFloat,2} = phi.eval,
    T::Type{<:AbstractFloat} = Float64,
    threads = false,
)
    n_equations = maximum_globalIndex(phi)

    x = zeros(T, n_equations)

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            for j in 1:mesh.m1
                if phi.onoff[i,j]
                    id = phi.gIndex[i,j]
                    x[id] = phisolution[i,j]
                end
            end
        end
    elseif !threads
        for i in 1:mesh.l1
            for j in 1:mesh.m1
                if phi.onoff[i,j]
                    id = phi.gIndex[i,j]
                    x[id] = phisolution[i,j]
                end
            end
        end
    end

    return x
end

function vector_to_phi!(
    phi::CSPhi3D,
    mesh::UnionCSMesh3D;
    phisolution::Array{<:AbstractFloat,3} = phi.eval,
    T::Type{<:AbstractFloat} = Float64,
    threads = false,
)
    n_equations = maximum_globalIndex(phi)

    x = zeros(T, n_equations)

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            for j in 1:mesh.m1
                for k in 1:mesh.n1
                    if phi.onoff[i,j,k]
                        id = phi.gIndex[i,j,k]
                        x[id] = phisolution[i,j,k]
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
                        x[id] = phisolution[i,j,k]
                    end
                end
            end
        end
    end

    return x
end
