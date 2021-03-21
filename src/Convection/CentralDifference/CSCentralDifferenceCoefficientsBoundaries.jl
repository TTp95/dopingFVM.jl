"""

"""
function _convection_centralDifference_bounds_ end

function _convection_centralDifference_bounds_(
    i::Signed,
    velocityU::Array{AbstractFloat,1},
    phi::CSPhi1D,
    phiBounds::BoundsStructured,
    mesh::UnionCSMesh1D,
)
    ab = 0.0
    bb = 0.0

    if (phiBounds.cord == 'w') || (phiBounds.cord == 'W')
        velf = abs(velocityU[i])
        if (velocityU[i] >= 0.0)
            ab = 0.0
            bb = -1.0 * velf * phiBounds.ρ * (1.0) * phiBounds.eval
        elseif (velocityU[i] <= 0.0)
            ab = velf * phiBounds.ρ * (1.0)
            bb = 0.0
        else (velocityU[i] == 0.0)
            ab = 0.0
            bb = 0.0
        end
    elseif (phiBounds.cord == 'e') || (phiBounds.cord == 'E')
        velf = abs(velocityU[i+1])
        if (velocityU[i+1] <= 0.0)
            ab = 0.0
            bb = -1.0 * velf * phiBounds.ρ * (1.0) * phiBounds.eval
        elseif (velocityU[i+1] >= 0.0)
            ab = velf * phiBounds.ρ * (1.0)
            bb = 0.0
        else (velocityU[i+1] == 0.0)
            ab = 0.0
            bb = 0.0
        end
    else
        error("Bound error in coordinate definition: $(phiBounds.cord)")
    end

    return ab, bb
end

function _diffusion_centralDifference_bounds_(
    i::Signed,
    j::Signed,
    velocityU::Array{AbstractFloat,2},
    velocityV::Array{AbstractFloat,2},
    phi::CSPhi2D,
    phiBounds::BoundsStructured,
    mesh::UnionCSMesh2D,
)
    ab = 0.0
    bb = 0.0

    if (phiBounds.cord == 'w') || (phiBounds.cord == 'W')
        velf = abs(velocityU[i,j])
        if (velocityU[i,j] >= 0.0)
            ab = 0.0
            bb = -1.0 * velf * phiBounds.ρ * (mesh.dy[j]) * phiBounds.eval
        elseif (velocityU[i,j] <= 0.0)
            ab = velf * phiBounds.ρ * (mesh.dy[j])
            bb = 0.0
        else (velocityU[i,j] == 0.0)
            ab = 0.0
            bb = 0.0
        end
    elseif (phiBounds.cord == 'e') || (phiBounds.cord == 'E')
        velf = abs(velocityU[i+1,j])
        if (velocityU[i+1,j] <= 0.0)
            ab = 0.0
            bb = -1.0 * velf * phiBounds.ρ * (mesh.dy[j]) * phiBounds.eval
        elseif (velocityU[i+1,j] >= 0.0)
            ab = phiBounds.ρ * (mesh.dy[j])
            bb = 0.0
        else (velocityU[i+1,j] == 0.0)
            ab = 0.0
            bb = 0.0
        end
    elseif (phiBounds.cord == 's') || (phiBounds.cord == 'S')
        velf = abs(velocityV[i,j-1])
        if (velocityV[i,j-1] >= 0.0)
            ab = 0.0
            bb = -1.0 * velf * phiBounds.ρ * (mesh.dx[i]) * phiBounds.eval
        elseif (velocityV[i,j-1] <= 0.0)
            ab = velf * phiBounds.ρ * (mesh.dx[i])
            bb = 0.0
        else (velocityV[i,j-1] == 0.0)
            ab = 0.0
            bb = 0.0
        end
    elseif (phiBounds.cord == 'n') || (phiBounds.cord == 'N')
        velf = abs(velocityV[i,j+1])
        if (velocityV[i,j+1] <= 0.0)
            ab = 0.0
            bb = -1.0 * velf * phiBounds.ρ * (mesh.dx[i]) * phiBounds.eval
        elseif (velocityV[i,j+1] >= 0.0)
            ab = velf * phiBounds.ρ * (mesh.dx[i])
            bb = 0.0
        else (velocityV[i,j+1] == 0.0)
            ab = 0.0
            bb = 0.0
        end
    else
        error("Bound error in coordinate definition: $(phiBounds.cord)")
    end

    return ab, bb
end

function _diffusion_centralDifference_bounds_(
    i::Signed,
    j::Signed,
    k::Signed,
    velocityU::Array{AbstractFloat,3},
    velocityV::Array{AbstractFloat,3},
    velocityW::Array{AbstractFloat,3},
    phi::CSPhi3D,
    phiBounds::BoundsStructured,
    mesh::UnionCSMesh3D,
)
    ab = 0.0
    bb = 0.0

    if (phiBounds.cord == 'w') || (phiBounds.cord == 'W')
        velf = abs(velocityU[i,j,k])
        if (velocityU[i,j,k] >= 0.0)
            ab = 0.0
            bb = -1.0 * velf * phiBounds.ρ * (mesh.dy[j] * mesh.dz[k]) * phiBounds.eval
        elseif (velocityU[i,j,k] <= 0.0)
            ab = velf * phiBounds.ρ * (mesh.dy[j] * mesh.dz[k])
            bb = 0.0
        else (velocityU[i,j,k] == 0.0)
            ab = 0.0
            bb = 0.0
        end
    elseif (phiBounds.cord == 'e') || (phiBounds.cord == 'E')
        velf = abs(velocityU[i+1,j,k])
        if (velocityU[i+1,j,k] <= 0.0)
            ab = 0.0
            bb = -1.0 * velf * phiBounds.ρ * (mesh.dy[j] * mesh.dz[k]) * phiBounds.eval
        elseif (velocityU[i+1,j,k] >= 0.0)
            ab = phiBounds.ρ * (mesh.dy[j] * mesh.dz[k])
            bb = 0.0
        else (velocityU[i+1,j,k] == 0.0)
            ab = 0.0
            bb = 0.0
        end
    elseif (phiBounds.cord == 's') || (phiBounds.cord == 'S')
        velf = abs(velocityV[i,j-1,k])
        if (velocityV[i,j-1,k] >= 0.0)
            ab = 0.0
            bb = -1.0 * velf * phiBounds.ρ * (mesh.dx[i] * mesh.dz[k]) * phiBounds.eval
        elseif (velocityV[i,j-1,k] <= 0.0)
            ab = velf * phiBounds.ρ * (mesh.dx[i] * mesh.dz[k])
            bb = 0.0
        else (velocityV[i,j-1,k] == 0.0)
            ab = 0.0
            bb = 0.0
        end
    elseif (phiBounds.cord == 'n') || (phiBounds.cord == 'N')
        velf = abs(velocityV[i,j+1,k])
        if (velocityV[i,j+1,k] <= 0.0)
            ab = 0.0
            bb = -1.0 * velf * phiBounds.ρ * (mesh.dx[i] * mesh.dz[k]) * phiBounds.eval
        elseif (velocityV[i,j+1,k] >= 0.0)
            ab = velf * phiBounds.ρ * (mesh.dx[i] * mesh.dz[k])
            bb = 0.0
        else (velocityV[i,j+1,k] == 0.0)
            ab = 0.0
            bb = 0.0
        end
    elseif (phiBounds.cord == 'b') || (phiBounds.cord == 'B')
        velf = abs(velocityV[i,j,k-1])
        if (velocityV[i,j,k-1] >= 0.0)
            ab = 0.0
            bb = -1.0 * velf * phiBounds.ρ * (mesh.dx[i] * mesh.dy[j]) * phiBounds.eval
        elseif (velocityV[i,j,k-1] <= 0.0)
            ab = velf * phiBounds.ρ * (mesh.dx[i] * mesh.dy[j])
            bb = 0.0
        else (velocityV[i,j,k-1] == 0.0)
            ab = 0.0
            bb = 0.0
        end
    elseif (phiBounds.cord == 't') || (phiBounds.cord == 'T')
        velf = abs(velocityV[i,j,k+1])
        if (velocityV[i,j,k+1] <= 0.0)
            ab = 0.0
            bb = -1.0 * velf * phiBounds.ρ * (mesh.dx[i] * mesh.dy[j]) * phiBounds.eval
        elseif (velocityV[i,j,k+1] >= 0.0)
            ab = velf * phiBounds.ρ * (mesh.dx[i] * mesh.dy[j])
            bb = 0.0
        else (velocityV[i,j,k+1] == 0.0)
            ab = 0.0
            bb = 0.0
        end
    else
        error("Bound error in coordinate definition: $(phiBounds.cord)")
    end

    return ab, bb
end
