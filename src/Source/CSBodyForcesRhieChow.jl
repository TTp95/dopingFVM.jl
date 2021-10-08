"""

"""
function discretize_bodyForcesRhieChow end

function discretize_bodyForcesRhieChow(
    velocity::CSVelocity1D,
    mesh::UnionCSMesh1D;
    T::Type{<:AbstractFloat} = Float64,
)
    n_equations = maximum_globalIndex(velocity.u)

    # b = spzeros(T, n_equations)
    b = zeros(T, n_equations)

    for i in 1:mesh.l1
        if velocity.u.onoff[i]
            id = velocity.u.gIndex[i]

            if velocity.u.onoff[i-1] && velocity.u.onoff[i+1]
                d0 = 0.5 * mesh.dx[i-1]
                d1 = 0.5 * mesh.dx[i]
                d2 = 0.5 * mesh.dx[i+1]

                gf1 = d0 / (d0 + d1)
                num = velocity.u.sourceC[i-1] * d1 + velocity.u.sourceC[i] * d0
                den =  d0 + d1
                a1 = (num / den) * (d0 + d1)

                gf2 = d2 / (d1 + d2)
                num = velocity.u.sourceC[i] * d2 + velocity.u.sourceC[i+1] * d1
                den =  d1 + d2
                a2 = (num / den) * (d1 + d2)

                @inbounds b[id]  = (a2 - a1) * 1.0

            else
                @inbounds b[id]  = velocity.u.sourceC[i] * mesh.vol[i]
            end

        end
    end

    return b
end

function discretize_bodyForcesRhieChow(
    velocity::CSVelocity2D,
    mesh::UnionCSMesh2D;
    T::Type{<:AbstractFloat} = Float64,
)
    nx_equations = maximum_globalIndex(velocity.u)
    ny_equations = maximum_globalIndex(velocity.v)

    # bx = spzeros(T, nx_equations)
    # by = spzeros(T, ny_equations)
    bx = zeros(T, nx_equations)
    by = zeros(T, ny_equations)

    for i in 1:mesh.l1
        for j in 1:mesh.m1
            id = velocity.u.gIndex[i,j]

            if velocity.u.onoff[i,j]
                if velocity.u.onoff[i-1,j] && velocity.u.onoff[i+1,j]
                    d0 = 0.5 * mesh.dx[i-1]
                    d1 = 0.5 * mesh.dx[i]
                    d2 = 0.5 * mesh.dx[i+1]

                    gf1 = d0 / (d0 + d1)
                    num = velocity.u.sourceC[i-1,j] * d1 + velocity.u.sourceC[i,j] * d0
                    den =  d0 + d1
                    a1 = (num / den) * (d0 + d1)

                    gf2 = d2 / (d1 + d2)
                    num = velocity.u.sourceC[i,j] * d2 + velocity.u.sourceC[i+1,j] * d1
                    den =  d1 + d2
                    a2 = (num / den) * (d1 + d2)

                    @inbounds bx[id]  = (a2 - a1) * mesh.dy[j]

                else
                    @inbounds bx[id]  = velocity.u.sourceC[i,j] * mesh.vol[i,j]
                end
            end

            if velocity.v.onoff[i,j]
                if velocity.v.onoff[i,j-1] && velocity.v.onoff[i,j+1]
                    d0 = 0.5 * mesh.dy[j-1]
                    d1 = 0.5 * mesh.dy[j]
                    d2 = 0.5 * mesh.dy[j+1]

                    gf1 = d0 / (d0 + d1)
                    num = velocity.v.sourceC[i,j-1] * d1 + velocity.v.sourceC[i,j] * d0
                    den =  d0 + d1
                    a1 = (num / den) * (d0 + d1)

                    gf2 = d2 / (d1 + d2)
                    num = velocity.v.sourceC[i,j] * d2 + velocity.v.sourceC[i,j+1] * d1
                    den =  d1 + d2
                    a2 = (num / den) * (d1 + d2)

                    @inbounds by[id]  = (a2 - a1) * mesh.dx[i]

                else
                    @inbounds by[id]  = velocity.v.sourceC[i,j] * mesh.vol[i,j]
                end
            end

        end
    end

    return A, b
end

function discretize_bodyForcesRhieChow(
    velocity::CSVelocity3D,
    mesh::UnionCSMesh3D;
    T::Type{<:AbstractFloat} = Float64,
)
    nx_equations = maximum_globalIndex(velocity.u)
    ny_equations = maximum_globalIndex(velocity.v)
    nz_equations = maximum_globalIndex(velocity.w)

    # bx = spzeros(T, nx_equations)
    # by = spzeros(T, ny_equations)
    # bz = spzeros(T, nz_equations)
    bx = zeros(T, nx_equations)
    by = zeros(T, ny_equations)
    bz = zeros(T, nz_equations)

    for i in 1:mesh.l1
        for j in 1:mesh.m1
            for k in 1:mesh.n1

                if velocity.u.onoff[i,j,k]
                    id = velocity.u.gIndex[i,j,k]

                    if velocity.u.onoff[i-1,j,k] && velocity.u.onoff[i+1,j,k]
                        d0 = 0.5 * mesh.dx[i-1]
                        d1 = 0.5 * mesh.dx[i]
                        d2 = 0.5 * mesh.dx[i+1]

                        gf1 = d0 / (d0 + d1)
                        num = velocity.u.sourceC[i-1,j,k] * d1 + velocity.u.sourceC[i,j,k] * d0
                        den =  d0 + d1
                        a1 = (num / den) * (d0 + d1)

                        gf2 = d2 / (d1 + d2)
                        num = velocity.u.sourceC[i,j,k] * d2 + velocity.u.sourceC[i+1,j,k] * d1
                        den =  d1 + d2
                        a2 = (num / den) * (d1 + d2)

                        @inbounds bx[id]  = (a2 - a1) * (mesh.dy[j] * mesh.dz[k])

                    else
                        @inbounds bx[id]  = velocity.u.sourceC[i,j,k] * mesh.vol[i,j,k]
                    end
                end

                if velocity.v.onoff[i,j,k]
                    if velocity.v.onoff[i,j-1,k] && velocity.v.onoff[i,j+1,k]
                        d0 = 0.5 * mesh.dy[j-1]
                        d1 = 0.5 * mesh.dy[j]
                        d2 = 0.5 * mesh.dy[j+1]

                        gf1 = d0 / (d0 + d1)
                        num = velocity.v.sourceC[i,j-1,k] * d1 + velocity.v.sourceC[i,j,k] * d0
                        den =  d0 + d1
                        a1 = (num / den) * (d0 + d1)

                        gf2 = d2 / (d1 + d2)
                        num = velocity.v.sourceC[i,j,k] * d2 + velocity.v.sourceC[i,j+1,k] * d1
                        den =  d1 + d2
                        a2 = (num / den) * (d1 + d2)

                        @inbounds by[id]  = (a2 - a1) * (mesh.dx[i] * mesh.dz[k])

                    else
                        @inbounds by[id]  = velocity.v.sourceC[i,j,k] * mesh.vol[i,j,k]
                    end
                end

                if velocity.w.onoff[i,j,k]
                    if velocity.w.onoff[i,j,k-1] && velocity.w.onoff[i,j,k+1]
                        d0 = 0.5 * mesh.dz[k-1]
                        d1 = 0.5 * mesh.dz[k]
                        d2 = 0.5 * mesh.dz[k+1]

                        gf1 = d0 / (d0 + d1)
                        num = velocity.w.sourceC[i,j,k-1] * d1 + velocity.w.sourceC[i,j,k] * d0
                        den =  d0 + d1
                        a1 = (num / den) * (d0 + d1)

                        gf2 = d2 / (d1 + d2)
                        num = velocity.w.sourceC[i,j,k] * d2 + velocity.w.sourceC[i,j,k+1] * d1
                        den =  d1 + d2
                        a2 = (num / den) * (d1 + d2)

                        @inbounds bz[id]  = (a2 - a1) * (mesh.dx[i] * mesh.dy[j])

                    else
                        @inbounds bz[id]  = velocity.w.sourceC[i,j,k] * mesh.vol[i,j,k]
                    end
                end

            end
        end
    end

    return A, b
end
