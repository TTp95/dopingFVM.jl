"""

"""
function array_gradient end

function array_gradient(
    arrayG::Array{<:AbstractFloat,1},
    phi::CSPhi1D,
    mesh::UnionCSMesh1D;
    T::Type{<:AbstractFloat} = Float64,
    threads = false,
)
    n_equations = maximum_globalIndex(phi)

    bx = zeros(T, n_equations)

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
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
    elseif !threads
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
    end

    return bx
end

function array_gradient(
    arrayG::Array{<:AbstractFloat,2},
    phi::CSPhi2D,
    mesh::UnionCSMesh2D;
    T::Type{<:AbstractFloat} = Float64,
    threads = false,
)
    n_equations = maximum_globalIndex(phi)

    bx = zeros(T, n_equations)
    by = zeros(T, n_equations)

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            for j in 1:mesh.m1
                if phi.onoff[i,j]
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
    elseif !threads
        for i in 1:mesh.l1
            for j in 1:mesh.m1
                if phi.onoff[i,j]
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
    threads = false,
)
    n_equations = maximum_globalIndex(phi)

    bx = zeros(T, n_equations)
    by = zeros(T, n_equations)
    bz = zeros(T, n_equations)

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            for j in 1:mesh.m1
                for k in 1:mesh.n1
                    if phi.onoff[i,j,k]
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
    elseif !threads
        for i in 1:mesh.l1
            for j in 1:mesh.m1
                for k in 1:mesh.n1
                    if phi.onoff[i,j,k]
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
