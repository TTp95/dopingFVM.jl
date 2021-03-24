"""

"""
function vector_to_phi! end

function vector_to_phi!(
    x::Vector{<:AbstractFloat},
    phi::CSPhi1D,
    mesh::UnionCSMesh1D;
    phisolution::Array{<:AbstractFloat,1} = phi.eval,
    T::Type{<:AbstractFloat} = Float64,
    threads = false,
)

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            if phi.onoff[i]
                id = phi.gIndex[i]
                phisolution[i] = x[id]
            end
        end
    elseif !threads
        for i in 1:mesh.l1
            if phi.onoff[i]
                id = phi.gIndex[i]
                phisolution[i] = x[id]
            end
        end
    end

    return nothing
end

function vector_to_phi!(
    x::Vector{<:AbstractFloat},
    phi::CSPhi2D,
    mesh::UnionCSMesh2D;
    phisolution::Array{<:AbstractFloat,2} = phi.eval,
    T::Type{<:AbstractFloat} = Float64,
    threads = false,
)

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            for j in 1:mesh.m1
                if phi.onoff[i,j]
                    id = phi.gIndex[i,j]
                    phisolution[i,j] = x[id]
                end
            end
        end
    elseif !threads
        for i in 1:mesh.l1
            for j in 1:mesh.m1
                if phi.onoff[i,j]
                    id = phi.gIndex[i,j]
                    phisolution[i,j] = x[id]
                end
            end
        end
    end

    return nothing
end

function vector_to_phi!(
    x::Vector{<:AbstractFloat},
    phi::CSPhi3D,
    mesh::UnionCSMesh3D;
    phisolution::Array{<:AbstractFloat,3} = phi.eval,
    T::Type{<:AbstractFloat} = Float64,
    threads = false,
)

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            for j in 1:mesh.m1
                for k in 1:mesh.n1
                    if phi.onoff[i,j,k]
                        id = phi.gIndex[i,j,k]
                        phisolution[i,j,k] = x[id]
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
                        phisolution[i,j,k] = x[id]
                    end
                end
            end
        end
    end

    return nothing
end
