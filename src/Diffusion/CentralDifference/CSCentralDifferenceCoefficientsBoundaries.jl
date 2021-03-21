"""

"""
function _diffusion_centralDifference_bounds_ end

function _diffusion_centralDifference_bounds_(
    i::Signed,
    phi::CSPhi1D,
    phiBounds::BoundsStructured,
    mesh::UnionCSMesh1D,
)
    ab = 0.0
    bb = 0.0

    d1 = 0.5 * mesh.dx[i]

    if (phiBounds.kind == 1) || (phiBounds.kind == 5) || (phiBounds.kind == 6)|| (phiBounds.kind == 7)
        ab = phiBounds.Γ * (1.0 / d1)
        bb = ab * phiBounds.a

    elseif (phiBounds.kind == 2)
        ab = 0.0
        bb = phiBounds.a * 1.0

    elseif (phiBounds.kind == 3)
        aux = phiBounds.Γ / d1
        num = phiBounds.a * (aux) * 1.0
        den = phiBounds.a + aux
        Req = num / den
        ab = Req
        bb = Req * phiBounds.b

    elseif (phiBounds.kind == 4) || (phiBounds.kind == 8)
        ab = 0.0
        bb = 0.0

    else
        error("Boundary condition id = $(phiBounds.kind) unimplemented")

    end

    return ab, bb
end

function _diffusion_centralDifference_bounds_(
    i::Signed,
    j::Signed,
    phi::CSPhi2D,
    phiBounds::BoundsStructured,
    mesh::UnionCSMesh2D,
)
    ab = 0.0
    bb = 0.0

    if (phiBounds.cord == 'w') || (phiBounds.cord == 'e') || (phiBounds.cord == 'W') || (phiBounds.cord == 'E')
        d1 = 0.5 * mesh.dx[i]
        area = mesh.dy[j]
    elseif (phiBounds.cord == 's') || (phiBounds.cord == 'n') || (phiBounds.cord == 'S') || (phiBounds.cord == 'N')
        d1 = 0.5 * mesh.dy[j]
        area = mesh.dx[i]
    else #Default case
        d1 = 0.5 * mesh.dx[i]
        area = mesh.dy[j]
    end

    if (phiBounds.kind == 1) || (phiBounds.kind == 5) || (phiBounds.kind == 6)|| (phiBounds.kind == 7)
        ab = phiBounds.Γ * (area / d1)
        bb = ab * phiBounds.a

    elseif (phiBounds.kind == 2)
        ab = 0.0
        bb = phiBounds.a * area

    elseif (phiBounds.kind == 3)
        aux = phiBounds.Γ / d1
        num = phiBounds.a * (aux) * area
        den = phiBounds.a + aux
        Req = (num / den)
        ab = Req
        bb = Req * phiBounds.b

    elseif (phiBounds.kind == 4) || (phiBounds.kind == 8)
        ab = 0.0
        bb = 0.0

    else
        error("Boundary condition id = $(phiBounds.kind) unimplemented")

    end

    return ab, bb
end

function _diffusion_centralDifference_bounds_(
    i::Signed,
    j::Signed,
    k::Signed,
    phi::CSPhi3D,
    phiBounds::BoundsStructured,
    mesh::UnionCSMesh3D,
)
    ab = 0.0
    bb = 0.0

    if (phiBounds.cord == 'w') || (phiBounds.cord == 'e') || (phiBounds.cord == 'W') || (phiBounds.cord == 'E')
        d1 = 0.5 * mesh.dx[i]
        area = mesh.dy[j] * mesh.dz[k]
    elseif (phiBounds.cord == 's') || (phiBounds.cord == 'n') || (phiBounds.cord == 'S') || (phiBounds.cord == 'N')
        d1 = 0.5 * mesh.dy[j]
        area = mesh.dx[i] * mesh.dz[k]
    elseif (phiBounds.cord == 'b') || (phiBounds.cord == 't') || (phiBounds.cord == 'B') || (phiBounds.cord == 'T')
        d1 = 0.5 * mesh.dz[k]
        area = mesh.dx[i] * mesh.dy[j]
    else #Default case
        d1 = 0.5 * mesh.dx[i]
        area = mesh.dy[j] * mesh.dz[k]
    end

    if (phiBounds.kind == 1) || (phiBounds.kind == 5) || (phiBounds.kind == 6)|| (phiBounds.kind == 7)
        ab = phiBounds.Γ * (area / d1)
        bb = ab * phiBounds.a

    elseif (phiBounds.kind == 2)
        ab = 0.0
        bb = phiBounds.a * area

    elseif (phiBounds.kind == 3)
        aux = phiBounds.Γ / d1
        num = phiBounds.a * (aux) * area
        den = phiBounds.a + aux
        Req = num / den
        ab = Req
        bb = Req * phiBounds.b

    elseif (phiBounds.kind == 4) || (phiBounds.kind == 8)
        ab = 0.0
        bb = 0.0

    else
        error("Boundary condition id = $(phiBounds.kind) unimplemented")

    end

    return ab, bb
end
