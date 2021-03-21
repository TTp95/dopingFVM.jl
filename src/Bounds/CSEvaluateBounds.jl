"""

"""
function _evaluate_bounds!_ end

function _evaluate_bounds!_(
    i::Signed,
    phi::CSPhi1D,
    mesh::UnionCSMesh1D,
    Boundsvalue::BoundsStructured;
)
    if (kind == 1) || (kind == 5) || (kind == 6) || (kind == 7)
        Boundsvalue.eval = Boundsvalue.a

    elseif (kind == 2)
        Diff = 1.0 / (0.5 * mesh.dx[i])
        num = Boundsvalue.Γ * Diff * phi.eval[i] - Boundsvalue.a
        den = Boundsvalue.Γ * Diff
        Boundsvalue.eval = num / den

    elseif (kind == 3)
        Diff = Boundsvalue.a / mesh.dx[i]
        num = (Boundsvalue.a * Boundsvalue.b) + (Diff * phi.eval[i])
        den = Boundsvalue.a + Diff
        Boundsvalue.eval = num / den


    elseif (kind == 4) || (kind == 1)

        if (mesh.l1 == 1) && ((Boundsvalue.cord == 'w') || (Boundsvalue.cord == 'W'))
            Boundsvalue.eval = phi.eval[i]

        elseif ((Boundsvalue.cord == 'w') || (Boundsvalue.cord == 'W')) && (i == mesh.l1)
            Boundsvalue.eval = phi.eval[i]

        elseif ((Boundsvalue.cord == 'w') || (Boundsvalue.cord == 'W')) && !phi.onoff[i+1]
            Boundsvalue.eval = phi.eval[i]

        elseif ((Boundsvalue.cord == 'w') || (Boundsvalue.cord == 'W')) && phi.onoff[i+1]
            num = phi.eval[i+1] - phi.eval[i]
            den = mesh.x[i+1] - mesh.x[i]
            grad =  a / b
            Boundsvalue.eval = phi.eval[i] - (grad * 0.5 * mesh.dx[i])

        end

        if (mesh.l1 == 1) && ((Boundsvalue.cord == 'e') || (Boundsvalue.cord == 'E'))
            Boundsvalue.eval = phi.eval[i]

        elseif ((Boundsvalue.cord == 'e') || (Boundsvalue.cord == 'E')) && (i == 1)
            Boundsvalue.eval = phi.eval[i]

        elseif ((Boundsvalue.cord == 'e') || (Boundsvalue.cord == 'E')) && !phi.onoff[i-1]
            Boundsvalue.eval = phi.eval[i]

        elseif ((Boundsvalue.cord == 'e') || (Boundsvalue.cord == 'E')) && phi.onoff[i-1]
            num = phi.eval[i] - phi.eval[i-1]
            den = mesh.x[i] - mesh.x[i-1]
            grad =  a / b
            Boundsvalue.eval = phi.eval[i] + (grad * 0.5 * mesh.dx[i])

        end

    else
        error("Boundary condition $(Boundsvalue.kind) unimplemented...")
    end

    return nothing
end

function _evaluate_bounds!_(
    i::Signed,
    j::Signed,
    phi::CSPhi2D,
    mesh::UnionCSMesh2D,
    Boundsvalue::BoundsStructured;
)
    if (phiBounds.cord == 'w') || (phiBounds.cord == 'e') || (phiBounds.cord == 'W') || (phiBounds.cord == 'E')
        dh = mesh.dx[i]
    elseif (phiBounds.cord == 's') || (phiBounds.cord == 'n') || (phiBounds.cord == 'S') || (phiBounds.cord == 'N')
        dh = mesh.dy[j]
    else #Default case
        error("Boundary of node i = $(i), j = $(j) - error in cord: $(Boundsvalue.cord)")
    end

    if (kind == 1) || (kind == 5) || (kind == 6) || (kind == 7)
        Boundsvalue.eval = Boundsvalue.a

    elseif (kind == 2)
        Diff = 1.0 / (0.5 * dh)
        num = Boundsvalue.Γ * Diff * phi.eval[i,j] - Boundsvalue.a
        den = Boundsvalue.Γ * Diff
        Boundsvalue.eval = num / den

    elseif (kind == 3)
        Diff = Boundsvalue.a / dh
        num = (Boundsvalue.a * Boundsvalue.b) + (Diff * phi.eval[i,j])
        den = Boundsvalue.a + Diff
        Boundsvalue.eval = num / den


    elseif (kind == 4) || (kind == 1)

        if (mesh.l1 == 1) && ((Boundsvalue.cord == 'w') || (Boundsvalue.cord == 'W'))
            Boundsvalue.eval = phi.eval[i,j]

        elseif ((Boundsvalue.cord == 'w') || (Boundsvalue.cord == 'W')) && (i == mesh.l1)
            Boundsvalue.eval = phi.eval[i,j]

        elseif ((Boundsvalue.cord == 'w') || (Boundsvalue.cord == 'W')) && !phi.onoff[i+1,j]
            Boundsvalue.eval = phi.eval[i,j]

        elseif ((Boundsvalue.cord == 'w') || (Boundsvalue.cord == 'W')) && phi.onoff[i+1,j]
            num = phi.eval[i+1,j] - phi.eval[i,j]
            den = mesh.x[i+1] - mesh.x[i]
            grad =  a / b
            Boundsvalue.eval = phi.eval[i,j] - (grad * 0.5 * mesh.dx[i])

        end


        if (mesh.l1 == 1) && ((Boundsvalue.cord == 'e') || (Boundsvalue.cord == 'E'))
            Boundsvalue.eval = phi.eval[i,j]

        elseif ((Boundsvalue.cord == 'e') || (Boundsvalue.cord == 'E')) && (i == 1)
            Boundsvalue.eval = phi.eval[i,j]

        elseif ((Boundsvalue.cord == 'e') || (Boundsvalue.cord == 'E')) && !phi.onoff[i-1,j]
            Boundsvalue.eval = phi.eval[i,j]

        elseif ((Boundsvalue.cord == 'e') || (Boundsvalue.cord == 'E')) && phi.onoff[i-1,j]
            num = phi.eval[i,j] - phi.eval[i-1,j]
            den = mesh.x[i] - mesh.x[i-1]
            grad =  a / b
            Boundsvalue.eval = phi.eval[i,j] + (grad * 0.5 * mesh.dx[i])

        end


        if (mesh.m1 == 1) && ((Boundsvalue.cord == 's') || (Boundsvalue.cord == 'S'))
            Boundsvalue.eval = phi.eval[i,j]

        elseif ((Boundsvalue.cord == 's') || (Boundsvalue.cord == 'S')) && (j == mesh.m1)
            Boundsvalue.eval = phi.eval[i,j]

        elseif ((Boundsvalue.cord == 's') || (Boundsvalue.cord == 'S')) && !phi.onoff[i,j+1]
            Boundsvalue.eval = phi.eval[i,j]

        elseif ((Boundsvalue.cord == 's') || (Boundsvalue.cord == 'S')) && phi.onoff[i,j+1]
            num = phi.eval[i,j+1] - phi.eval[i,j]
            den = mesh.y[j+1] - mesh.y[j]
            grad =  a / b
            Boundsvalue.eval = phi.eval[i,j] - (grad * 0.5 * mesh.dy[j])

        end


        if (mesh.m1 == 1) && ((Boundsvalue.cord == 'n') || (Boundsvalue.cord == 'N'))
            Boundsvalue.eval = phi.eval[i,j]

        elseif ((Boundsvalue.cord == 'n') || (Boundsvalue.cord == 'N')) && (j == 1)
            Boundsvalue.eval = phi.eval[i,j]

        elseif ((Boundsvalue.cord == 'n') || (Boundsvalue.cord == 'N')) && !phi.onoff[i,j-1]
            Boundsvalue.eval = phi.eval[i,j]

        elseif ((Boundsvalue.cord == 'n') || (Boundsvalue.cord == 'N')) && phi.onoff[i,j-1]
            num = phi.eval[i,j] - phi.eval[i,j-1]
            den = mesh.y[j] - mesh.y[j-1]
            grad =  a / b
            Boundsvalue.eval = phi.eval[i,j] + (grad * 0.5 * mesh.dy[j])

        end


    else
        error("Boundary condition $(Boundsvalue.kind) unimplemented...")
    end

    return nothing
end

function _evaluate_bounds!_(
    i::Signed,
    j::Signed,
    k::Signed,
    phi::CSPhi3D,
    mesh::UnionCSMesh3D,
    Boundsvalue::BoundsStructured;
)
    if (phiBounds.cord == 'w') || (phiBounds.cord == 'e') || (phiBounds.cord == 'W') || (phiBounds.cord == 'E')
        dh = mesh.dx[i]
    elseif (phiBounds.cord == 's') || (phiBounds.cord == 'n') || (phiBounds.cord == 'S') || (phiBounds.cord == 'N')
        dh = mesh.dy[j]
    elseif (phiBounds.cord == 'b') || (phiBounds.cord == 'B') || (phiBounds.cord == 't') || (phiBounds.cord == 'T')
        dh = mesh.dz[k]
    else #Default case
        error("Boundary of node i = $(i), j = $(j) - error in cord: $(Boundsvalue.cord)")
    end

    if (kind == 1) || (kind == 5) || (kind == 6) || (kind == 7)
        Boundsvalue.eval = Boundsvalue.a

    elseif (kind == 2)
        Diff = 1.0 / (0.5 * dh)
        num = Boundsvalue.Γ * Diff * phi.eval[i,j,k] - Boundsvalue.a
        den = Boundsvalue.Γ * Diff
        Boundsvalue.eval = num / den

    elseif (kind == 3)
        Diff = Boundsvalue.a / dh
        num = (Boundsvalue.a * Boundsvalue.b) + (Diff * phi.eval[i,j,k])
        den = Boundsvalue.a + Diff
        Boundsvalue.eval = num / den


    elseif (kind == 4) || (kind == 1)

        if (mesh.l1 == 1) && ((Boundsvalue.cord == 'w') || (Boundsvalue.cord == 'W'))
            Boundsvalue.eval = phi.eval[i,j,k]

        elseif ((Boundsvalue.cord == 'w') || (Boundsvalue.cord == 'W')) && (i == mesh.l1)
            Boundsvalue.eval = phi.eval[i,j,k]

        elseif ((Boundsvalue.cord == 'w') || (Boundsvalue.cord == 'W')) && !phi.onoff[i+1,j,k]
            Boundsvalue.eval = phi.eval[i,j,k]

        elseif ((Boundsvalue.cord == 'w') || (Boundsvalue.cord == 'W')) && phi.onoff[i+1,j,k]
            num = phi.eval[i+1,j,k] - phi.eval[i,j,k]
            den = mesh.x[i+1] - mesh.x[i]
            grad =  a / b
            Boundsvalue.eval = phi.eval[i,j,k] - (grad * 0.5 * mesh.dx[i])

        end


        if (mesh.l1 == 1) && ((Boundsvalue.cord == 'e') || (Boundsvalue.cord == 'E'))
            Boundsvalue.eval = phi.eval[i,j,k]

        elseif ((Boundsvalue.cord == 'e') || (Boundsvalue.cord == 'E')) && (i == 1)
            Boundsvalue.eval = phi.eval[i,j,k]

        elseif ((Boundsvalue.cord == 'e') || (Boundsvalue.cord == 'E')) && !phi.onoff[i-1,j,k]
            Boundsvalue.eval = phi.eval[i,j,k]

        elseif ((Boundsvalue.cord == 'e') || (Boundsvalue.cord == 'E')) && phi.onoff[i-1,j,k]
            num = phi.eval[i,j,k] - phi.eval[i-1,j,k]
            den = mesh.x[i] - mesh.x[i-1]
            grad =  a / b
            Boundsvalue.eval = phi.eval[i,j,k] + (grad * 0.5 * mesh.dx[i])

        end


        if (mesh.m1 == 1) && ((Boundsvalue.cord == 's') || (Boundsvalue.cord == 'S'))
            Boundsvalue.eval = phi.eval[i,j,k]

        elseif ((Boundsvalue.cord == 's') || (Boundsvalue.cord == 'S')) && (j == mesh.m1)
            Boundsvalue.eval = phi.eval[i,j,k]

        elseif ((Boundsvalue.cord == 's') || (Boundsvalue.cord == 'S')) && !phi.onoff[i,j+1,k]
            Boundsvalue.eval = phi.eval[i,j,k]

        elseif ((Boundsvalue.cord == 's') || (Boundsvalue.cord == 'S')) && phi.onoff[i,j+1,k]
            num = phi.eval[i,j+1,k] - phi.eval[i,j,k]
            den = mesh.y[j+1] - mesh.y[j]
            grad =  a / b
            Boundsvalue.eval = phi.eval[i,j,k] - (grad * 0.5 * mesh.dy[j])

        end


        if (mesh.m1 == 1) && ((Boundsvalue.cord == 'n') || (Boundsvalue.cord == 'N'))
            Boundsvalue.eval = phi.eval[i,j,k]

        elseif ((Boundsvalue.cord == 'n') || (Boundsvalue.cord == 'N')) && (j == 1)
            Boundsvalue.eval = phi.eval[i,j,k]

        elseif ((Boundsvalue.cord == 'n') || (Boundsvalue.cord == 'N')) && !phi.onoff[i,j-1,k]
            Boundsvalue.eval = phi.eval[i,j,k]

        elseif ((Boundsvalue.cord == 'n') || (Boundsvalue.cord == 'N')) && phi.onoff[i,j-1,k]
            num = phi.eval[i,j,k] - phi.eval[i,j-1,k]
            den = mesh.y[j] - mesh.y[j-1]
            grad =  a / b
            Boundsvalue.eval = phi.eval[i,j,k] + (grad * 0.5 * mesh.dy[j])

        end


        if (mesh.n1 == 1) && ((Boundsvalue.cord == 'b') || (Boundsvalue.cord == 'B'))
            Boundsvalue.eval = phi.eval[i,j,k]

        elseif ((Boundsvalue.cord == 'b') || (Boundsvalue.cord == 'B')) && (k == mesh.n1)
            Boundsvalue.eval = phi.eval[i,j,k]

        elseif ((Boundsvalue.cord == 'b') || (Boundsvalue.cord == 'B')) && !phi.onoff[i,j,k+1]
            Boundsvalue.eval = phi.eval[i,j,k]

        elseif ((Boundsvalue.cord == 'b') || (Boundsvalue.cord == 'B')) && phi.onoff[i,j,k+1]
            num = phi.eval[i,j,k+1] - phi.eval[i,j,k]
            den = mesh.z[k+1] - mesh.z[k]
            grad =  a / b
            Boundsvalue.eval = phi.eval[i,j,k] - (grad * 0.5 * mesh.dz[k])

        end


        if (mesh.n1 == 1) && ((Boundsvalue.cord == 'n') || (Boundsvalue.cord == 'N'))
            Boundsvalue.eval = phi.eval[i,j,k]

        elseif ((Boundsvalue.cord == 'n') || (Boundsvalue.cord == 'N')) && (k == 1)
            Boundsvalue.eval = phi.eval[i,j,k]

        elseif ((Boundsvalue.cord == 'n') || (Boundsvalue.cord == 'N')) && !phi.onoff[i,j,k-1]
            Boundsvalue.eval = phi.eval[i,j,k]

        elseif ((Boundsvalue.cord == 'n') || (Boundsvalue.cord == 'N')) && phi.onoff[i,j,k-1]
            num = phi.eval[i,j,k] - phi.eval[i,j,k-1]
            den = mesh.z[k] - mesh.z[k-1]
            grad =  a / b
            Boundsvalue.eval = phi.eval[i,j,k] + (grad * 0.5 * mesh.dz[k])

        end


    else
        error("Boundary condition $(Boundsvalue.kind) unimplemented...")
    end

    return nothing
end
