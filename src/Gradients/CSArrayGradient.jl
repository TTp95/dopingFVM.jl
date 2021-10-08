"""

"""
function array_gradient end

function array_gradient(
    arrayG::Array{<:AbstractFloat,1},
    phi::CSPhi1D,
    mesh::UnionCSMesh1D;
    T::Type{<:AbstractFloat} = Float64,
)
    n_equations = maximum_globalIndex(phi)

    bx = zeros(T, n_equations)

    for i in 1:mesh.l1
        if phi.onoff[i]
            if (i == 1)
                id = phi.gIndex[i]
                num = arrayG[i+1] - arrayG[i]
                den = 0.5 * (mesh.dx[i] + mesh.dx[i+1])
                bx[id] = (num / den) * mesh.vol[i]
            elseif (i == mesh.l1)
                id = phi.gIndex[i]
                num = arrayG[i] - arrayG[i-1]
                den = 0.5 * (mesh.dx[i] + mesh.dx[i-1])
                bx[id] = (num / den) * mesh.vol[i]
            elseif (phi.onoff[i-1]) && (phi.onoff[i+1])
                id = phi.gIndex[i]
                dx0 = 0.5 * mesh.dx[i-1]
                dx1 = 0.5 * mesh.dx[i]
                dx2 = 0.5 * mesh.dx[i+1]
                num = arrayG[i-1] * dx1 + arrayG[i] * dx0
                den = dx1 + dx0
                a = num / den
                num = arrayG[i] * dx2 + arrayG[i+1] * dx1
                den = dx1 + dx2
                b = num / den
                bx[id] = ((b - a) / mesh.dx[i]) * mesh.vol[i]
            elseif (!phi.onoff[i-1])
                id = phi.gIndex[i]
                num = arrayG[i+1] - arrayG[i]
                den = 0.5 * (mesh.dx[i] + mesh.dx[i+1])
                bx[id] = (num / den) * mesh.vol[i]
            elseif (!phi.onoff[i+1])
                id = phi.gIndex[i]
                num = arrayG[i] - arrayG[i-1]
                den = 0.5 * (mesh.dx[i] + mesh.dx[i-1])
                bx[id] = (num / den) * mesh.vol[i]
            else
                error("CONFIGURATION ERROR")
            end
        end
    end

    return bx
end

function array_gradient(
    arrayG::Array{<:AbstractFloat,2},
    phi::CSPhi2D,
    mesh::UnionCSMesh2D;
    T::Type{<:AbstractFloat} = Float64,
)
    n_equations = maximum_globalIndex(phi)

    bx = zeros(T, n_equations)
    by = zeros(T, n_equations)

    for i in 1:mesh.l1
        for j in 1:mesh.m1
            if phi.onoff[i,j]
                if (mesh.l1 != 1)
                    if (i == 1)
                        id = phi.gIndex[i,j]
                        num = arrayG[i+1,j] - arrayG[i,j]
                        den = 0.5 * (mesh.dx[i] + mesh.dx[i+1])
                        bx[id] = (num / den) * mesh.vol[i,j]
                    elseif (i == mesh.l1)
                        id = phi.gIndex[i,j]
                        num = arrayG[i,j] - arrayG[i-1,j]
                        den = 0.5 * (mesh.dx[i] + mesh.dx[i-1])
                        bx[id] = (num / den) * mesh.vol[i,j]
                    elseif (phi.onoff[i-1,j]) && (phi.onoff[i+1,j])
                        id = phi.gIndex[i,j]
                        dx0 = 0.5 * mesh.dx[i-1]
                        dx1 = 0.5 * mesh.dx[i]
                        dx2 = 0.5 * mesh.dx[i+1]
                        num = arrayG[i-1,j] * dx1 + arrayG[i,j] * dx0
                        den = dx1 + dx0
                        a = num / den
                        num = arrayG[i,j] * dx2 + arrayG[i+1,j] * dx1
                        den = dx1 + dx2
                        b = num / den
                        bx[id] = ((b - a) / mesh.dx[i]) * mesh.vol[i,j]
                    elseif (!phi.onoff[i-1,j])
                        id = phi.gIndex[i,j]
                        num = arrayG[i+1,j] - arrayG[i,j]
                        den = 0.5 * (mesh.dx[i] + mesh.dx[i+1])
                        bx[id] = (num / den) * mesh.vol[i,j]
                    elseif (!phi.onoff[i+1,j])
                        id = phi.gIndex[i,j]
                        num = arrayG[i,j] - arrayG[i-1,j]
                        den = 0.5 * (mesh.dx[i] + mesh.dx[i-1])
                        bx[id] = (num / den) * mesh.vol[i,j]
                    else
                        error("CONFIGURATION ERROR")
                    end
                end

                if (mesh.m1 != 1)
                    if (j == 1)
                        id = phi.gIndex[i,j]
                        num = arrayG[i,j+1] - arrayG[i,j]
                        den = 0.5 * (mesh.dy[j] + mesh.dy[j+1])
                        by[id] = (num / den) * mesh.vol[i,j]
                    elseif (j == mesh.m1)
                        id = phi.gIndex[i,j]
                        num = arrayG[i,j] - arrayG[i,j-1]
                        den = 0.5 * (mesh.dy[j] + mesh.dy[j-1])
                        by[id] = (num / den) * mesh.vol[i,j]
                    elseif (phi.onoff[i,j-1]) && (phi.onoff[i,j+1])
                        id = phi.gIndex[i,j]
                        dy0 = 0.5 * mesh.dy[j-1]
                        dy1 = 0.5 * mesh.dy[j]
                        dy2 = 0.5 * mesh.dy[j+1]
                        num = arrayG[i,j-1] * dy1 + arrayG[i,j] * dy0
                        den = dy1 + dy0
                        a = num / den
                        num = arrayG[i,j] * dy2 + arrayG[i,j+1] * dy1
                        den = dy1 + dy2
                        b = num / den
                        by[id] = ((b - a) / mesh.dy[j]) * mesh.vol[i,j]
                    elseif (!phi.onoff[i,j-1])
                        id = phi.gIndex[i,j]
                        num = arrayG[i,j+1] - arrayG[i,j]
                        den = 0.5 * (mesh.dy[j] + mesh.dy[j+1])
                        by[id] = (num / den) * mesh.vol[i,j]
                    elseif (!phi.onoff[i,j+1])
                        id = phi.gIndex[i,j]
                        num = arrayG[i,j] - arrayG[i,j-1]
                        den = 0.5 * (mesh.dy[j] + mesh.dy[j-1])
                        by[id] = (num / den) * mesh.vol[i,j]
                    else
                        error("CONFIGURATION ERROR")
                    end
                end
            end
        end
    end

    return bx, by
end

function array_gradient(
    arrayG::Array{<:AbstractFloat,3},
    phi::CSPhi3D,
    mesh::UnionCSMesh3D;
    T::Type{<:AbstractFloat} = Float64,
)
    n_equations = maximum_globalIndex(phi)

    bx = zeros(T, n_equations)
    by = zeros(T, n_equations)
    bz = zeros(T, n_equations)

    for i in 1:mesh.l1
        for j in 1:mesh.m1
            for k in 1:mesh.n1
                if phi.onoff[i,j,k]
                    if (mesh.l1 != 1)
                        if (i == 1)
                            id = phi.gIndex[i,j,k]
                            num = arrayG[i+1,j,k] - arrayG[i,j,k]
                            den = 0.5 * (mesh.dx[i] + mesh.dx[i+1])
                            bx[id] = (num / den) * mesh.vol[i,j,k]
                        elseif (i == mesh.l1)
                            id = phi.gIndex[i,j,k]
                            num = arrayG[i,j,k] - arrayG[i-1,j,k]
                            den = 0.5 * (mesh.dx[i] + mesh.dx[i-1])
                            bx[id] = (num / den) * mesh.vol[i,j,k]
                        elseif (phi.onoff[i-1,j,k]) && (phi.onoff[i+1,j,k])
                            id = phi.gIndex[i,j,k]
                            dx0 = 0.5 * mesh.dx[i-1]
                            dx1 = 0.5 * mesh.dx[i]
                            dx2 = 0.5 * mesh.dx[i+1]
                            num = arrayG[i-1,j,k] * dx1 + arrayG[i,j,k] * dx0
                            den = dx1 + dx0
                            a = num / den
                            num = arrayG[i,j,k] * dx2 + arrayG[i+1,j,k] * dx1
                            den = dx1 + dx2
                            b = num / den
                            bx[id] = ((b - a) / mesh.dx[i]) * mesh.vol[i,j,k]
                        elseif (!phi.onoff[i-1,j,k])
                            id = phi.gIndex[i,j,k]
                            num = arrayG[i+1,j,k] - arrayG[i,j,k]
                            den = 0.5 * (mesh.dx[i] + mesh.dx[i+1])
                            bx[id] = (num / den) * mesh.vol[i,j,k]
                        elseif (!phi.onoff[i+1,j,k])
                            id = phi.gIndex[i,j,k]
                            num = arrayG[i,j,k] - arrayG[i-1,j,k]
                            den = 0.5 * (mesh.dx[i] + mesh.dx[i-1])
                            bx[id] = (num / den) * mesh.vol[i,j,k]
                        else
                            error("CONFIGURATION ERROR")
                        end
                    end

                    if (mesh.m1 != 1)
                        if (j == 1)
                            id = phi.gIndex[i,j,k]
                            num = arrayG[i,j+1,k] - arrayG[i,j,k]
                            den = 0.5 * (mesh.dy[j] + mesh.dy[j+1])
                            by[id] = (num / den) * mesh.vol[i,j,k]
                        elseif (j == mesh.m1)
                            id = phi.gIndex[i,j,k]
                            num = arrayG[i,j,k] - arrayG[i,j-1,k]
                            den = 0.5 * (mesh.dy[j] + mesh.dy[j-1])
                            by[id] = (num / den) * mesh.vol[i,j,k]
                        elseif (phi.onoff[i,j-1,k]) && (phi.onoff[i,j+1,k])
                            id = phi.gIndex[i,j,k]
                            dy0 = 0.5 * mesh.dy[j-1]
                            dy1 = 0.5 * mesh.dy[j]
                            dy2 = 0.5 * mesh.dy[j+1]
                            num = arrayG[i,j-1,k] * dy1 + arrayG[i,j,k] * dy0
                            den = dy1 + dy0
                            a = num / den
                            num = arrayG[i,j,k] * dy2 + arrayG[i,j+1,k] * dy1
                            den = dy1 + dy2
                            b = num / den
                            by[id] = ((b - a) / mesh.dy[j]) * mesh.vol[i,j,k]
                        elseif (!phi.onoff[i,j-1,k])
                            id = phi.gIndex[i,j,k]
                            num = arrayG[i,j+1,k] - arrayG[i,j,k]
                            den = 0.5 * (mesh.dy[j] + mesh.dy[j+1])
                            by[id] = (num / den) * mesh.vol[i,j,k]
                        elseif (!phi.onoff[i,j+1,k])
                            id = phi.gIndex[i,j,k]
                            num = arrayG[i,j,k] - arrayG[i,j-1,k]
                            den = 0.5 * (mesh.dy[j] + mesh.dy[j-1])
                            by[id] = (num / den) * mesh.vol[i,j,k]
                        else
                            error("CONFIGURATION ERROR")
                        end
                    end

                    if (mesh.n1 != 1)
                        if (k == 1)
                            id = phi.gIndex[i,j,k]
                            num = arrayG[i,j,k+1] - arrayG[i,j,k]
                            den = 0.5 * (mesh.dz[k] + mesh.dz[k+1])
                            bz[id] = (num / den) * mesh.vol[i,j,k]
                        elseif (k == mesh.n1)
                            id = phi.gIndex[i,j,k]
                            num = arrayG[i,j,k] - arrayG[i,j,k-1]
                            den = 0.5 * (mesh.dz[k] + mesh.dz[k-1])
                            bz[id] = (num / den) * mesh.vol[i,j,k]
                        elseif (phi.onoff[i,j,k-1]) && (phi.onoff[i,j,k+1])
                            id = phi.gIndex[i,j,k]
                            dz0 = 0.5 * mesh.dz[k-1]
                            dz1 = 0.5 * mesh.dz[k]
                            dz2 = 0.5 * mesh.dz[k+1]
                            num = arrayG[i,j,k-1] * dz1 + arrayG[i,j,k] * dz0
                            den = dz1 + dz0
                            a = num / den
                            num = arrayG[i,j,k] * dz2 + arrayG[i,j,k+1] * dz1
                            den = dz1 + dz2
                            b = num / den
                            bz[id] = ((b - a) / mesh.dz[k]) * mesh.vol[i,j,k]
                        elseif (!phi.onoff[i,j,k-1])
                            id = phi.gIndex[i,j,k]
                            num = arrayG[i,j,k+1] - arrayG[i,j,k]
                            den = 0.5 * (mesh.dz[k] + mesh.dz[k+1])
                            bz[id] = (num / den) * mesh.vol[i,j,k]
                        elseif (!phi.onoff[i,j,k+1])
                            id = phi.gIndex[i,j,k]
                            num = arrayG[i,j,k] - arrayG[i,j,k-1]
                            den = 0.5 * (mesh.dz[k] + mesh.dz[k-1])
                            bz[id] = (num / den) * mesh.vol[i,j,k]
                        else
                            error("CONFIGURATION ERROR")
                        end
                    end

                end
            end
        end
    end

    return bx, by, bz
end

function array_gradient(
    arrayG::Array{<:AbstractFloat,1},
    phi::CSPhi1D,
    bounds::Dict{String,BoundsStructured},
    mesh::UnionCSMesh1D;
    T::Type{<:AbstractFloat} = Float64,
)
    n_equations = maximum_globalIndex(phi)

    bx = zeros(T, n_equations)

    for i in 1:mesh.l1
        if phi.onoff[i]
            id = phi.gIndex[i]

            #Auxiliar variables
            a = 0.0
            b = 0.0
            num = 0.0
            den = 0.0
            dx0 = 0.0
            dx1 = 0.0
            dx2 = 0.0
            boundvalue = 0.0

            if (i == 1)
                if phi.bounds[i]
                    dx1 = 0.5 * mesh.dx[i]
                    dx2 = 0.5 * mesh.dx[i+1]
                    num = arrayG[i] * dx2 + arrayG[i+1] * dx1
                    den = dx1 + dx2
                    b = num / den
                    boundvalue = find_bondValue(i, id, phi, mesh, bounds, 'w')
                    num = b - boundvalue
                    den = mesh.dx[i]
                    bx[id] = (num / den) * mesh.vol[i]
                else
                    num = arrayG[i+1] - arrayG[i]
                    den = 0.5 * (mesh.dx[i+1] + mesh.dx[i])
                    bx[id] = (num / den) * mesh.vol[i]
                end
            elseif (i == mesh.l1)
                if phi.bounds[i]
                    dx0 = 0.5 * mesh.dx[i-1]
                    dx1 = 0.5 * mesh.dx[i]
                    num = arrayG[i] * dx0 + arrayG[i-1] * dx1
                    den = dx0 + dx1
                    a = num / den
                    boundvalue = find_bondValue(i, id, phi, mesh, bounds, 'e')
                    num = boundvalue - a
                    den = mesh.dx[i]
                    bx[id] = (num / den) * mesh.vol[i]
                else
                    num = arrayG[i] - arrayG[i-1]
                    den = 0.5 * (mesh.dx[i] + mesh.dx[i-1])
                    bx[id] = (num / den) * mesh.vol[i]
                end
            elseif (phi.onoff[i-1]) && (phi.onoff[i+1])
                dx0 = 0.5 * mesh.dx[i-1]
                dx1 = 0.5 * mesh.dx[i]
                dx2 = 0.5 * mesh.dx[i+1]
                num = arrayG[i-1] * dx1 + arrayG[i] * dx0
                den = dx1 + dx0
                a = num / den
                num = arrayG[i] * dx2 + arrayG[i+1] * dx1
                den = dx1 + dx2
                b = num / den
                bx[id] = ((b - a) / mesh.dx[i]) * mesh.vol[i]
            elseif (!phi.onoff[i-1])
                dx1 = 0.5 * mesh.dx[i]
                dx2 = 0.5 * mesh.dx[i+1]
                num = arrayG[i] * dx2 + arrayG[i+1] * dx1
                den = dx1 + dx2
                b = num / den
                if phi.bounds[i]
                    boundvalue = find_bondValue(i, id, phi, mesh, bounds, 'w')
                    num = b - boundvalue
                    den = mesh.dx[i]
                    bx[id] = (num / den) * mesh.vol[i]
                else
                    num = b - arrayG[i-1]
                    den =mesh.dx[i]
                    bx[id] = (num / den) * mesh.vol[i]
                end
            elseif (!phi.onoff[i+1])
                dx0 = 0.5 * mesh.dx[i-1]
                dx1 = 0.5 * mesh.dx[i]
                num = arrayG[i-1] * dx1 + arrayG[i] * dx0
                den = dx1 + dx0
                a = num / den
                if phi.bounds[i]
                    boundvalue = find_bondValue(i, id, phi, mesh, bounds, 'e')
                    num = boundvalue - a
                    den = mesh.dx[i]
                    bx[id] = (num / den) * mesh.vol[i]
                else
                    num = arrayG[i+1] - a
                    den = mesh.dx[i]
                    bx[id] = (num / den) * mesh.vol[i]
                end
            else
                error("CONFIGURATION ERROR")
            end
        end
    end

    return bx
end

function array_gradient(
    arrayG::Array{<:AbstractFloat,2},
    phi::CSPhi2D,
    bounds::Dict{String,BoundsStructured},
    mesh::UnionCSMesh2D;
    T::Type{<:AbstractFloat} = Float64,
)
    n_equations = maximum_globalIndex(phi)

    bx = zeros(T, n_equations)
    by = zeros(T, n_equations)

    for i in 1:mesh.l1
        for j in 1:mesh.m1
            if phi.onoff[i,j]
                id = phi.gIndex[i,j]

                #Auxiliar variables
                a = 0.0
                b = 0.0
                num = 0.0
                den = 0.0
                dx0 = 0.0
                dx1 = 0.0
                dx2 = 0.0
                dy0 = 0.0
                dy1 = 0.0
                dy2 = 0.0
                boundvalue = 0.0

                if (mesh.l1 != 1)
                    if (i == 1)
                        if phi.bounds[i,j]
                            dx1 = 0.5 * mesh.dx[i]
                            dx2 = 0.5 * mesh.dx[i+1]
                            num = arrayG[i,j] * dx2 + arrayG[i+1,j] * dx1
                            den = dx1 + dx2
                            b = num / den
                            boundvalue = find_bondValue(i, j, id, phi, mesh, bounds, 'w')
                            num = b - boundvalue
                            den = mesh.dx[i]
                            bx[id] = (num / den) * mesh.vol[i,j]
                        else
                            num = arrayG[i+1,j] - arrayG[i,j]
                            den = 0.5 * (mesh.dx[i+1] + mesh.dx[i])
                            bx[id] = (num / den) * mesh.vol[i,j]
                        end
                    elseif (i == mesh.l1)
                        if phi.bounds[i,j]
                            dx0 = 0.5 * mesh.dx[i-1]
                            dx1 = 0.5 * mesh.dx[i]
                            num = arrayG[i,j] * dx0 + arrayG[i-1,j] * dx1
                            den = dx0 + dx1
                            a = num / den
                            boundvalue = find_bondValue(i, j, id, phi, mesh, bounds, 'e')
                            num = boundvalue - a
                            den = mesh.dx[i]
                            bx[id] = (num / den) * mesh.vol[i,j]
                        else
                            num = arrayG[i,j] - arrayG[i-1,j]
                            den = 0.5 * (mesh.dx[i] + mesh.dx[i-1])
                            bx[id] = (num / den) * mesh.vol[i,j]
                        end
                    elseif (phi.onoff[i-1,j]) && (phi.onoff[i+1,j])
                        dx0 = 0.5 * mesh.dx[i-1]
                        dx1 = 0.5 * mesh.dx[i]
                        dx2 = 0.5 * mesh.dx[i+1]
                        num = arrayG[i-1,j] * dx1 + arrayG[i,j] * dx0
                        den = dx1 + dx0
                        a = num / den
                        num = arrayG[i,j] * dx2 + arrayG[i+1,j] * dx1
                        den = dx1 + dx2
                        b = num / den
                        bx[id] = ((b - a) / mesh.dx[i]) * mesh.vol[i,j]
                    elseif (!phi.onoff[i-1,j])
                        dx1 = 0.5 * mesh.dx[i]
                        dx2 = 0.5 * mesh.dx[i+1]
                        num = arrayG[i,j] * dx2 + arrayG[i+1,j] * dx1
                        den = dx1 + dx2
                        b = num / den
                        if phi.bounds[i,j]
                            boundvalue = find_bondValue(i, j, id, phi, mesh, bounds, 'w')
                            num = b - boundvalue
                            den = mesh.dx[i]
                            bx[id] = (num / den) * mesh.vol[i,j]
                        else
                            num = b - arrayG[i-1,j]
                            den =mesh.dx[i]
                            bx[id] = (num / den) * mesh.vol[i,j]
                        end
                    elseif (!phi.onoff[i+1,j])
                        dx0 = 0.5 * mesh.dx[i-1]
                        dx1 = 0.5 * mesh.dx[i]
                        num = arrayG[i-1,j] * dx1 + arrayG[i,j] * dx0
                        den = dx1 + dx0
                        a = num / den
                        if phi.bounds[i,j]
                            boundvalue = find_bondValue(i, j, id, phi, mesh, bounds, 'e')
                            num = boundvalue - a
                            den = mesh.dx[i]
                            bx[id] = (num / den) * mesh.vol[i,j]
                        else
                            num = arrayG[i+1,j] - a
                            den = mesh.dx[i]
                            bx[id] = (num / den) * mesh.vol[i,j]
                        end
                    else
                        error("CONFIGURATION ERROR")
                    end
                end

                if (mesh.m1 != 1)
                    if (j == 1)
                        if phi.bounds[i,j]
                            dy1 = 0.5 * mesh.dy[j]
                            dy2 = 0.5 * mesh.dy[j+1]
                            num = arrayG[i,j] * dy2 + arrayG[i,j+1] * dy1
                            den = dy1 + dy2
                            b = num / den
                            boundvalue = find_bondValue(i, j, id, phi, mesh, bounds, 's')
                            num = b - boundvalue
                            den = mesh.dy[j]
                            by[id] = (num / den) * mesh.vol[i,j]
                        else
                            num = arrayG[i,j+1] - arrayG[i,j]
                            den = 0.5 * (mesh.dy[j+1] + mesh.dy[j])
                            by[id] = (num / den) * mesh.vol[i,j]
                        end
                    elseif (j == mesh.m1)
                        if phi.bounds[i,j]
                            dy0 = 0.5 * mesh.dy[j-1]
                            dy1 = 0.5 * mesh.dy[j]
                            num = arrayG[i,j] * dy0 + arrayG[i,j-1] * dy1
                            den = dy0 + dy1
                            a = num / den
                            boundvalue = find_bondValue(i, j, id, phi, mesh, bounds, 'n')
                            num = boundvalue - a
                            den = mesh.dy[j]
                            by[id] = (num / den) * mesh.vol[i,j]
                        else
                            num = arrayG[i,j] - arrayG[i,j-1]
                            den = 0.5 * (mesh.dy[j] + mesh.dy[j-1])
                            by[id] = (num / den) * mesh.vol[i,j]
                        end
                    elseif (phi.onoff[i,j-1]) && (phi.onoff[i,j+1])
                        dx0 = 0.5 * mesh.dy[j-1]
                        dx1 = 0.5 * mesh.dy[j]
                        dx2 = 0.5 * mesh.dy[j+1]
                        num = arrayG[i,j-1] * dy1 + arrayG[i,j] * dy0
                        den = dy1 + dy0
                        a = num / den
                        num = arrayG[i,j] * dy2 + arrayG[i,j+1] * dy1
                        den = dy1 + dy2
                        b = num / den
                        by[id] = ((b - a) / mesh.dy[j]) * mesh.vol[i,j]
                    elseif (!phi.onoff[i,j-1])
                        dy1 = 0.5 * mesh.dy[j]
                        dy2 = 0.5 * mesh.dy[j+1]
                        num = arrayG[i,j] * dy2 + arrayG[i,j+1] * dy1
                        den = dy1 + dy2
                        b = num / den
                        if phi.bounds[i,j]
                            boundvalue = find_bondValue(i, j, id, phi, mesh, bounds, 's')
                            num = b - boundvalue
                            den = mesh.dy[j]
                            by[id] = (num / den) * mesh.vol[i,j]
                        else
                            num = b - arrayG[i,j-1]
                            den = mesh.dy[j]
                            by[id] = (num / den) * mesh.vol[i,j]
                        end
                    elseif (!phi.onoff[i,j+1])
                        dy0 = 0.5 * mesh.dy[j-1]
                        dy1 = 0.5 * mesh.dy[j]
                        num = arrayG[i,j-1] * dy1 + arrayG[i,j] * dy0
                        den = dy1 + dy0
                        a = num / den
                        if phi.bounds[i,j]
                            boundvalue = find_bondValue(i, j, id, phi, mesh, bounds, 'n')
                            num = boundvalue - a
                            den = mesh.dy[j]
                            by[id] = (num / den) * mesh.vol[i,j]
                        else
                            num = arrayG[i,j+1] - a
                            den = mesh.dy[j]
                            by[id] = (num / den) * mesh.vol[i,j]
                        end
                    else
                        error("CONFIGURATION ERROR")
                    end
                end
            end
        end
    end

    return bx, by
end

function array_gradient(
    arrayG::Array{<:AbstractFloat,3},
    phi::CSPhi3D,
    bounds::Dict{String,BoundsStructured},
    mesh::UnionCSMesh3D;
    T::Type{<:AbstractFloat} = Float64,
)
    n_equations = maximum_globalIndex(phi)

    bx = zeros(T, n_equations)
    by = zeros(T, n_equations)
    bz = zeros(T, n_equations)

    for i in 1:mesh.l1
        for j in 1:mesh.m1
            for k in 1:mesh.n1
                if phi.onoff[i,j,k]
                    id = phi.gIndex[i,j,k]

                    #Auxiliar variables
                    a = 0.0
                    b = 0.0
                    num = 0.0
                    den = 0.0
                    dx0 = 0.0
                    dx1 = 0.0
                    dx2 = 0.0
                    dy0 = 0.0
                    dy1 = 0.0
                    dy2 = 0.0
                    dz0 = 0.0
                    dz1 = 0.0
                    dz2 = 0.0
                    boundvalue = 0.0

                    if (mesh.l1 != 1)
                        if (i == 1)
                            if phi.bounds[i,j,k]
                                dx1 = 0.5 * mesh.dx[i]
                                dx2 = 0.5 * mesh.dx[i+1]
                                num = arrayG[i,j,k] * dx2 + arrayG[i+1,j,k] * dx1
                                den = dx1 + dx2
                                b = num / den
                                boundvalue = find_bondValue(i, j, k, id, phi, mesh, bounds, 'w')
                                num = b - boundvalue
                                den = mesh.dx[i]
                                bx[id] = (num / den) * mesh.vol[i,j,k]
                            else
                                num = arrayG[i+1,j,k] - arrayG[i,j,k]
                                den = 0.5 * (mesh.dx[i+1] + mesh.dx[i])
                                bx[id] = (num / den) * mesh.vol[i,j,k]
                            end
                        elseif (i == mesh.l1)
                            if phi.bounds[i,j,k]
                                dx0 = 0.5 * mesh.dx[i-1]
                                dx1 = 0.5 * mesh.dx[i]
                                num = arrayG[i,j,k] * dx0 + arrayG[i-1,j,k] * dx1
                                den = dx0 + dx1
                                a = num / den
                                boundvalue = find_bondValue(i, j, k, id, phi, mesh, bounds, 'e')
                                num = boundvalue - a
                                den = mesh.dx[i]
                                bx[id] = (num / den) * mesh.vol[i,j,k]
                            else
                                num = arrayG[i,j,k] - arrayG[i-1,j,k]
                                den = 0.5 * (mesh.dx[i] + mesh.dx[i-1])
                                bx[id] = (num / den) * mesh.vol[i,j,k]
                            end
                        elseif (phi.onoff[i-1,j,k]) && (phi.onoff[i+1,j,k])
                            dx0 = 0.5 * mesh.dx[i-1]
                            dx1 = 0.5 * mesh.dx[i]
                            dx2 = 0.5 * mesh.dx[i+1]
                            num = arrayG[i-1,j,k] * dx1 + arrayG[i,j,k] * dx0
                            den = dx1 + dx0
                            a = num / den
                            num = arrayG[i,j,k] * dx2 + arrayG[i+1,j,k] * dx1
                            den = dx1 + dx2
                            b = num / den
                            bx[id] = ((b - a) / mesh.dx[i]) * mesh.vol[i,j,k]
                        elseif (!phi.onoff[i-1,j])
                            dx1 = 0.5 * mesh.dx[i]
                            dx2 = 0.5 * mesh.dx[i+1]
                            num = arrayG[i,j,k] * dx2 + arrayG[i+1,j,k] * dx1
                            den = dx1 + dx2
                            b = num / den
                            if phi.bounds[i,j,k]
                                boundvalue = find_bondValue(i, j, k, id, phi, mesh, bounds, 'w')
                                num = b - boundvalue
                                den = mesh.dx[i]
                                bx[id] = (num / den) * mesh.vol[i,j,k]
                            else
                                num = b - arrayG[i-1,j,k]
                                den = mesh.dx[i]
                                bx[id] = (num / den) * mesh.vol[i,j,k]
                            end
                        elseif (!phi.onoff[i+1,j,k])
                            dx0 = 0.5 * mesh.dx[i-1]
                            dx1 = 0.5 * mesh.dx[i]
                            num = arrayG[i-1,j,k] * dx1 + arrayG[i,j,k] * dx0
                            den = dx1 + dx0
                            a = num / den
                            if phi.bounds[i,j,k]
                                boundvalue = find_bondValue(i, j, k, id, phi, mesh, bounds, 'e')
                                num = boundvalue - a
                                den = mesh.dx[i]
                                bx[id] = (num / den) * mesh.vol[i,j,k]
                            else
                                num = arrayG[i+1,j,k] - a
                                den = mesh.dx[i]
                                bx[id] = (num / den) * mesh.vol[i,j,k]
                            end
                        else
                            error("CONFIGURATION ERROR")
                        end
                    end

                    if (mesh.m1 != 1)
                        if (j == 1)
                            if phi.bounds[i,j,k]
                                dy1 = 0.5 * mesh.dy[j]
                                dy2 = 0.5 * mesh.dy[j+1]
                                num = arrayG[i,j,k] * dy2 + arrayG[i,j+1,k] * dy1
                                den = dy1 + dy2
                                b = num / den
                                boundvalue = find_bondValue(i, j, k, id, phi, mesh, bounds, 's')
                                num = b - boundvalue
                                den = mesh.dy[j]
                                by[id] = (num / den) * mesh.vol[i,j,k]
                            else
                                num = arrayG[i,j+1,k] - arrayG[i,j,k]
                                den = 0.5 * (mesh.dy[j+1] + mesh.dy[j])
                                by[id] = (num / den) * mesh.vol[i,j,k]
                            end
                        elseif (j == mesh.m1)
                            if phi.bounds[i,j,k]
                                dy0 = 0.5 * mesh.dy[j-1]
                                dy1 = 0.5 * mesh.dy[j]
                                num = arrayG[i,j,k] * dy0 + arrayG[i,j-1,k] * dy1
                                den = dy0 + dy1
                                a = num / den
                                boundvalue = find_bondValue(i, j, k, id, phi, mesh, bounds, 'n')
                                num = boundvalue - a
                                den = mesh.dy[j]
                                by[id] = (num / den) * mesh.vol[i,j,k]
                            else
                                num = arrayG[i,j,k] - arrayG[i,j-1,k]
                                den = 0.5 * (mesh.dy[j] + mesh.dy[j-1])
                                by[id] = (num / den) * mesh.vol[i,j,k]
                            end
                        elseif (phi.onoff[i,j-1,k]) && (phi.onoff[i,j+1,k])
                            dx0 = 0.5 * mesh.dy[j-1]
                            dx1 = 0.5 * mesh.dy[j]
                            dx2 = 0.5 * mesh.dy[j+1]
                            num = arrayG[i,j-1,k] * dy1 + arrayG[i,j,k] * dy0
                            den = dy1 + dy0
                            a = num / den
                            num = arrayG[i,j,k] * dy2 + arrayG[i,j+1,k] * dy1
                            den = dy1 + dy2
                            b = num / den
                            by[id] = ((b - a) / mesh.dy[j]) * mesh.vol[i,j,k]
                        elseif (!phi.onoff[i,j-1,k])
                            dy1 = 0.5 * mesh.dy[j]
                            dy2 = 0.5 * mesh.dy[j+1]
                            num = arrayG[i,j,k] * dy2 + arrayG[i,j+1,k] * dy1
                            den = dy1 + dy2
                            b = num / den
                            if phi.bounds[i,j,k]
                                boundvalue = find_bondValue(i, j, k, id, phi, mesh, bounds, 's')
                                num = b - boundvalue
                                den = mesh.dy[j]
                                by[id] = (num / den) * mesh.vol[i,j,k]
                            else
                                num = b - arrayG[i,j-1,k]
                                den = mesh.dy[j]
                                by[id] = (num / den) * mesh.vol[i,j,k]
                            end
                        elseif (!phi.onoff[i,j+1,k])
                            dy0 = 0.5 * mesh.dy[j-1]
                            dy1 = 0.5 * mesh.dy[j]
                            num = arrayG[i,j-1,k] * dy1 + arrayG[i,j,k] * dy0
                            den = dy1 + dy0
                            a = num / den
                            if phi.bounds[i,j,k]
                                boundvalue = find_bondValue(i, j, k, id, phi, mesh, bounds, 'n')
                                num = boundvalue - a
                                den = mesh.dy[j]
                                by[id] = (num / den) * mesh.vol[i,j,k]
                            else
                                num = arrayG[i,j+1,k] - a
                                den = mesh.dy[j]
                                by[id] = (num / den) * mesh.vol[i,j,k]
                            end
                        else
                            error("CONFIGURATION ERROR")
                        end
                    end

                    if (mesh.n1 != 1)
                        if (k == 1)
                            if phi.bounds[i,j,k]
                                dz1 = 0.5 * mesh.dz[k]
                                dz2 = 0.5 * mesh.dz[k+1]
                                num = arrayG[i,j,k] * dz2 + arrayG[i,j,k+1] * dz1
                                den = dz1 + dz2
                                b = num / den
                                boundvalue = find_bondValue(i, j, k, id, phi, mesh, bounds, 'b')
                                num = b - boundvalue
                                den = mesh.dz[k]
                                bz[id] = (num / den) * mesh.vol[i,j,k]
                            else
                                num = arrayG[i,j,k+1] - arrayG[i,j,k]
                                den = 0.5 * (mesh.dz[k+1] + mesh.dz[k])
                                bz[id] = (num / den) * mesh.vol[i,j,k]
                            end
                        elseif (z == mesh.n1)
                            if phi.bounds[i,j,k]
                                dz0 = 0.5 * mesh.dz[k-1]
                                dz1 = 0.5 * mesh.dz[k]
                                num = arrayG[i,j,k] * dz0 + arrayG[i,j,k-1] * dz1
                                den = dz0 + dz1
                                a = num / den
                                boundvalue = find_bondValue(i, j, k, id, phi, mesh, bounds, 't')
                                num = boundvalue - a
                                den = mesh.dz[k]
                                bz[id] = (num / den) * mesh.vol[i,j,k]
                            else
                                num = arrayG[i,j,k] - arrayG[i,j,k-1]
                                den = 0.5 * (mesh.dz[k] + mesh.dz[k-1])
                                bz[id] = (num / den) * mesh.vol[i,j,k]
                            end
                        elseif (phi.onoff[i,j,k-1]) && (phi.onoff[i,j,k+1])
                            dz0 = 0.5 * mesh.dz[k-1]
                            dz1 = 0.5 * mesh.dz[k]
                            dz2 = 0.5 * mesh.dz[k+1]
                            num = arrayG[i,j,k-1] * dz1 + arrayG[i,j,k] * dz0
                            den = dz1 + dz0
                            a = num / den
                            num = arrayG[i,j,k] * dz2 + arrayG[i,j,k+1] * dz1
                            den = dz1 + dz2
                            b = num / den
                            bz[id] = ((b - a) / mesh.dz[k]) * mesh.vol[i,j,k]
                        elseif (!phi.onoff[i,j,k-1])
                            dz1 = 0.5 * mesh.dz[k]
                            dz2 = 0.5 * mesh.dz[k+1]
                            num = arrayG[i,j,k] * dz2 + arrayG[i,j,k+1] * dz1
                            den = dz1 + dz2
                            b = num / den
                            if phi.bounds[i,j,k]
                                boundvalue = find_bondValue(i, j, k, id, phi, mesh, bounds, 'b')
                                num = b - boundvalue
                                den = mesh.dz[k]
                                bz[id] = (num / den) * mesh.vol[i,j,k]
                            else
                                num = b - arrayG[i,j,k-1]
                                den = mesh.dz[k]
                                bz[id] = (num / den) * mesh.vol[i,j,k]
                            end
                        elseif (!phi.onoff[i,j,k+1])
                            dz0 = 0.5 * mesh.dz[k-1]
                            dz1 = 0.5 * mesh.dz[k]
                            num = arrayG[i,j,k-1] * dz1 + arrayG[i,j,k] * dz0
                            den = dz1 + dz0
                            a = num / den
                            if phi.bounds[i,j,k]
                                boundvalue = find_bondValue(i, j, k, id, phi, mesh, bounds, 't')
                                num = boundvalue - a
                                den = mesh.dz[k]
                                bz[id] = (num / den) * mesh.vol[i,j,k]
                            else
                                num = arrayG[i,j,k+1] - a
                                den = mesh.dz[k]
                                bz[id] = (num / den) * mesh.vol[i,j,k]
                            end
                        else
                            error("CONFIGURATION ERROR")
                        end
                    end
                end
            end
        end
    end

    return bx, by, bz
end
