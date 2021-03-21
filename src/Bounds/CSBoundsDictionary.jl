"""

"""
function create_BoundsDict end

function create_BoundsDict(
    phi::CSPhi1D,
    mesh::UnionCSMesh1D;
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
    threads = false,
)

    dictBounds = Dict{String,BoundsStructured}()

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            if phi.bounds[i]
                for nn in 1:phi.nbounds[i]
                    gid = phi.gIndex[i]
                    @inbounds dictBounds["$(gid)g$(nn)"] = create_BoundsStructured(
                    ; T = T, N = N,
                    )
                end
            end
        end
    elseif !threads
        for i in 1:mesh.l1
            if phi.bounds[i]
                for nn in 1:phi.nbounds[i]
                    gid = phi.gIndex[i]
                    @inbounds dictBounds["$(gid)g$(nn)"] = create_BoundsStructured(
                    ; T = T, N = N,
                    )
                end
            end
        end
    end

    return dictBounds
end

function create_BoundsDict(
    phi::CSPhi2D,
    mesh::UnionCSMesh2D;
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
    threads = false,
)

    dictBounds = Dict{String,BoundsStructured}()

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            for j in 1:mesh.m1
                if phi.bounds[i,j]
                    for nn in 1:phi.nbounds[i,j]
                        gid = phi.gIndex[i,j]
                        @inbounds dictBounds["$(gid)g$(nn)"] = create_BoundsStructured(
                        ; T = T, N = N,
                        )
                    end
                end
            end
        end

    elseif !threads
        for i in 1:mesh.l1
            for j in 1:mesh.m1
                if phi.bounds[i,j]
                    for nn in 1:phi.nbounds[i,j]
                        gid = phi.gIndex[i,j]
                        @inbounds dictBounds["$(gid)g$(nn)"] = create_BoundsStructured(
                        ; T = T, N = N,
                        )
                    end
                end
            end
        end
    end

    return dictBounds
end

function create_BoundsDict(
    phi::CSPhi3D,
    mesh::UnionCSMesh3D;
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
    threads = false,
)

    dictBounds = Dict{String,BoundsStructured}()

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            for j in 1:mesh.m1
                for k in 1:mesh.n1
                    if phi.bounds[i,j,k]
                        for nn in 1:phi.nbounds[i,j,k]
                            gid = phi.gIndex[i,j,k]
                            @inbounds dictBounds["$(gid)g$(nn)"] = create_BoundsStructured(
                            ; T = T, N = N,
                            )
                        end
                    end
                end
            end
        end

    elseif !threads
        for i in 1:mesh.l1
            for j in 1:mesh.m1
                for k in 1:mesh.n1
                    if phi.bounds[i,j,k]
                        for nn in 1:phi.nbounds[i,j,k]
                            gid = phi.gIndex[i,j,k]
                            @inbounds dictBounds["$(gid)g$(nn)"] = create_BoundsStructured(
                            ; T = T, N = N,
                            )
                        end
                    end
                end
            end
        end

    end

    return dictBounds
end
