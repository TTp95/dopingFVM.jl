"""

"""
function _globalIndex_ end

function _globalIndex_(
    mesh::UnionCSMesh1D,
    i::Signed,
)
    gindex =  i
    return gindex
end

function _globalIndex_(
    mesh::UnionCSMesh2D,
    i::Signed,
    j::Signed,
)
    gindex = i + (j - 1) * mesh.l1
    return gindex
end

function _globalIndex_(
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
function assign_globalIndex! end

function assign_globalIndex!(
    phi::CSPhi1D,
    mesh::UnionCSMesh1D;
)
    jumps = 0

    for i in 1:mesh.l1
        if phi.onoff[i]
            @inbounds phi.gIndex[i] = _globalIndex_(mesh,i) - jumps
        elseif !phi.onoff[i]
            jumps += 1
            @inbounds phi.gIndex[i] = -1
        end
    end

    return nothing
end

function assign_globalIndex!(
    phi::CSPhi2D,
    mesh::UnionCSMesh2D;
)
    jumps = 0

    for j in 1:mesh.m1
        for i in 1:mesh.l1
            if phi.onoff[i,j]
                @inbounds phi.gIndex[i,j] = _globalIndex_(mesh,i,j) - jumps
            elseif !phi.onoff[i,j]
                jumps += 1
                @inbounds phi.gIndex[i,j] = -1
            end
        end
    end

    return nothing
end

function assign_globalIndex!(
    phi::CSPhi3D,
    mesh::UnionCSMesh3D;
)
    jumps = 0

    for k in 1:mesh.n1
        for j in 1:mesh.m1
            for i in 1:mesh.l1
                if phi.onoff[i,j,k]
                    @inbounds phi.gIndex[i,j,k] = _globalIndex_(mesh,i,j,k) - jumps
                elseif !phi.onoff[i,j,k]
                    jumps += 1
                    @inbounds phi.gIndex[i,j,k] = -1
                end
            end
        end
    end

    return nothing
end

function assign_globalIndex!(
    vel::CSVelocity1D,
    mesh::UnionCSMesh1D;
)
    assign_globalIndex!(vel.u, mesh)
    assign_globalIndex!(vel.p, mesh)

    return nothing
end

function assign_globalIndex!(
    vel::CSVelocity2D,
    mesh::UnionCSMesh2D;
)
    assign_globalIndex!(vel.u, mesh)
    assign_globalIndex!(vel.v, mesh)
    assign_globalIndex!(vel.p, mesh)

    return nothing
end

function assign_globalIndex!(
    vel::CSVelocity3D,
    mesh::UnionCSMesh3D;
)
    assign_globalIndex!(vel.u, mesh)
    assign_globalIndex!(vel.v, mesh)
    assign_globalIndex!(vel.w, mesh)
    assign_globalIndex!(vel.p, mesh)

    return nothing
end

"""

"""
function maximum_globalIndex end

function maximum_globalIndex(phi::CSPhi1D)
    return maximum(phi.gIndex)
end

function maximum_globalIndex(phi::CSPhi2D)
    return maximum(phi.gIndex)
end

function maximum_globalIndex(phi::CSPhi3D)
    return maximum(phi.gIndex)
end
