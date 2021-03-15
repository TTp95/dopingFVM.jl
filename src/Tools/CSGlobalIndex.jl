"""

"""
function _globalIndexCS_ end

function _globalIndexCS_(
    mesh::UnionCSMesh1D,
    i::Signed,
)
    gindex =  i
    return gindex
end

function _globalIndexCS_(
    mesh::UnionCSMesh2D,
    i::Signed,
    j::Signed,
)
    gindex = i + (j - 1) * mesh.l1
    return gindex
end

function _globalIndexCS_(
    mesh::UnionCSMesh3D,
    i::Signed,
    j::Signed,
    k::Signed,
)
    gindex = i + (j - 1) * mesh.l1 + (k - 1) * mesh.m1 * mesh.l1
    return gindex
end

"""

"""
function assign_globalIndexCS! end

function assign_globalIndexCS!(
    phi::CSPhi1D,
    mesh::UnionCSMesh1D;
    threads = false,
)
    jumps = 0

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            if phi.onoff[i]
                @inbounds phi.gIndex[i] = _globalIndexCS_(mesh,i) - jumps
            elseif !phi.onoff[i]
                jumps += 1
                @inbounds phi.gIndex[i] = -1
            end
        end
    elseif !threads
        for i in 1:mesh.l1
            if phi.onoff[i]
                @inbounds phi.gIndex[i] = _globalIndexCS_(mesh,i) - jumps
            elseif !phi.onoff[i]
                jumps += 1
                @inbounds phi.gIndex[i] = -1
            end
        end
    end

    return nothing
end

function assign_globalIndexCS!(
    phi::CSPhi2D,
    mesh::UnionCSMesh2D;
    threads = false,
)
    jumps = 0

    if threads
        Base.Threads.@threads for j in 1:mesh.m1
            for i in 1:mesh.l1
                if phi.onoff[i,j]
                    @inbounds phi.gIndex[i,j] = _globalIndexCS_(mesh,i,j) - jumps
                elseif !phi.onoff[i,j]
                    jumps += 1
                    @inbounds phi.gIndex[i,j] = -1
                end
            end
        end
    elseif !threads
        for j in 1:mesh.m1
            for i in 1:mesh.l1
                if phi.onoff[i,j]
                    @inbounds phi.gIndex[i,j] = _globalIndexCS_(mesh,i,j) - jumps
                elseif !phi.onoff[i,j]
                    jumps += 1
                    @inbounds phi.gIndex[i,j] = -1
                end
            end
        end
    end

    return nothing
end

function assign_globalIndexCS!(
    phi::CSPhi3D,
    mesh::UnionCSMesh3D;
    threads = false,
)
    jumps = 0

    if threads
        Base.Threads.@threads for k in 1:mesh.n1
            for j in 1:mesh.m1
                for i in 1:mesh.l1
                    if phi.onoff[i,j,k]
                        @inbounds phi.gIndex[i,j,k] = _globalIndexCS_(mesh,i,j,k) - jumps
                    elseif !phi.onoff[i,j,k]
                        jumps += 1
                        @inbounds phi.gIndex[i,j,k] = -1
                    end
                end
            end
        end
    elseif !threads
        for k in 1:mesh.n1
            for j in 1:mesh.m1
                for i in 1:mesh.l1
                    if phi.onoff[i,j,k]
                        @inbounds phi.gIndex[i,j,k] = _globalIndexCS_(mesh,i,j,k) - jumps
                    elseif !phi.onoff[i,j,k]
                        jumps += 1
                        @inbounds phi.gIndex[i,j,k] = -1
                    end
                end
            end
        end
    end

    return nothing
end

"""

"""
function maximum_globalIndexCS end

function maximum_globalIndexCS(phi::CSPhi1D)
    return maximum(phi.gIndex)
end

function maximum_globalIndexCS(phi::CSPhi2D)
    return maximum(phi.gIndex)
end

function maximum_globalIndexCS(phi::CSPhi3D)
    return maximum(phi.gIndex)
end
