"""

"""
function _projection_PC_pressure_laplacian_neighbors_(
    lenghts1::AbstractFloat,
    lenghts2::AbstractFloat,
    area::AbstractFloat,
)
    dCF = 0.5 * (lenghts1 + lenghts2)
    aF = -1.0 * (area / dCF)

    return aF
end

"""

"""
function _projection_PC_pressure_laplacian_ end

function _projection_PC_pressure_laplacian_(
    phi::CSPhi1D,
    mesh::UnionCSMesh1D;
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
)
    n_equations = maximum_globalIndex(phi)

    n = 0
    AI = zeros(N, (3 * n_equations))
    AJ = zeros(N, (3 * n_equations))
    AV = zeros(T, (3 * n_equations))

    for i in 1:mesh.l1
        if phi.onoff[i]
            id = phi.gIndex[i]

            #Auxiliar variables
            ac = 0.0
            aw = 0.0
            ae = 0.0

            #West Coefficents
            if (i != 1) && (mesh.l1 != 1) && (phi.onoff[i-1])
                @inbounds aw = _projection_PC_pressure_laplacian_neighbors_(
                    mesh.dx[i], mesh.dx[i-1], (1.0)
                )
                n += 1
                AI[n] = id
                AJ[n] = phi.gIndex[i-1]
                AV[n] = aw
            end

            #East Coefficents
            if (i != mesh.l1) && (mesh.l1 != 1)  && (phi.onoff[i+1])
                @inbounds ae = _projection_PC_pressure_laplacian_neighbors_(
                    mesh.dx[i], mesh.dx[i+1], (1.0)
                )
                n += 1
                AI[n] = id
                AJ[n] = phi.gIndex[i+1]
                AV[n] = ae
            end

            #Center Coefficent
            n += 1
            AI[n] = id
            AJ[n] = id
            AV[n] = -1.0 * (aw + ae)
        end
    end

    A = sparse(AI[1:n], AJ[1:n], AV[1:n])

    return A
end

function _projection_PC_pressure_laplacian_(
    phi::CSPhi2D,
    mesh::UnionCSMesh2D;
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
)
    n_equations = maximum_globalIndex(phi)

    n = 0
    AI = zeros(N, (5 * n_equations))
    AJ = zeros(N, (5 * n_equations))
    AV = zeros(T, (5 * n_equations))

    for i in 1:mesh.l1
        for j in 1:mesh.m1
            if phi.onoff[i,j]
                id = phi.gIndex[i,j]

                #Auxiliar variables
                ac = 0.0
                aw = 0.0
                ae = 0.0
                as = 0.0
                an = 0.0

                #West Coefficents
                if (i != 1) && (mesh.l1 != 1) && (phi.onoff[i-1,j])
                    @inbounds aw = _projection_PC_pressure_laplacian_neighbors_(
                        mesh.dx[i], mesh.dx[i-1], (mesh.dy[j])
                    )
                    n += 1
                    AI[n] = id
                    AJ[n] = phi.gIndex[i-1,j]
                    AV[n] = aw
                end

                #East Coefficents
                if (i != mesh.l1) && (mesh.l1 != 1) && (phi.onoff[i+1,j])
                    @inbounds ae = _projection_PC_pressure_laplacian_neighbors_(
                        mesh.dx[i], mesh.dx[i+1], (mesh.dy[j])
                    )
                    n += 1
                    AI[n] = id
                    AJ[n] = phi.gIndex[i+1,j]
                    AV[n] = ae
                end

                #South Coefficents
                if (j != 1) && (mesh.m1 != 1) && (phi.onoff[i,j-1])
                    @inbounds as = _projection_PC_pressure_laplacian_neighbors_(
                        mesh.dy[j], mesh.dy[j-1], (mesh.dx[i])
                    )
                    n += 1
                    AI[n] = id
                    AJ[n] = phi.gIndex[i,j-1]
                    AV[n] = as
                end

                #North Coefficents
                if (j != mesh.m1) && (mesh.m1 != 1) && (phi.onoff[i,j+1])
                    @inbounds an = _projection_PC_pressure_laplacian_neighbors_(
                        mesh.dy[j], mesh.dy[j+1], (mesh.dx[i])
                    )
                    n += 1
                    AI[n] = id
                    AJ[n] = phi.gIndex[i,j+1]
                    AV[n] = an
                end

                #Center Coefficent
                n += 1
                AI[n] = id
                AJ[n] = id
                AV[n] = -1.0 * (aw + ae + as + an)
            end
        end
    end

    A = sparse(AI[1:n], AJ[1:n], AV[1:n])

    return A
end

function _projection_PC_pressure_laplacian_(
    phi::CSPhi3D,
    mesh::UnionCSMesh3D;
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
)
    n_equations = maximum_globalIndex(phi)

    n = 0
    AI = zeros(N, (7 * n_equations))
    AJ = zeros(N, (7 * n_equations))
    AV = zeros(T, (7 * n_equations))

    for i in 1:mesh.l1
        for j in 1:mesh.m1
            for k in 1:mesh.n1
                if phi.onoff[i,j,k]
                    id = phi.gIndex[i,j,k]

                    #Auxiliar variables
                    ac = 0.0
                    aw = 0.0
                    ae = 0.0
                    as = 0.0
                    an = 0.0
                    ab = 0.0
                    at = 0.0

                    #West Coefficents
                    if (i != 1) && (mesh.l1 != 1) && (phi.onoff[i-1,j,k])
                        @inbounds aw = _projection_PC_pressure_laplacian_neighbors_(
                            mesh.dx[i], mesh.dx[i-1], (mesh.dy[j] * mesh.dz[k])
                        )
                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i-1,j,k]
                        AV[n] = aw
                    end

                    #East Coefficents
                    if (i != mesh.l1) && (mesh.l1 != 1) && (phi.onoff[i+1,j,k])
                        @inbounds ae = _projection_PC_pressure_laplacian_neighbors_(
                            mesh.dx[i], mesh.dx[i+1], (mesh.dy[j] * mesh.dz[k])
                        )
                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i+1,j,k]
                        AV[n] = ae
                    end

                    #South Coefficents
                    if (j != 1) && (mesh.m1 != 1) && (phi.onoff[i,j-1,k])
                        @inbounds as = _projection_PC_pressure_laplacian_neighbors_(
                            mesh.dy[j], mesh.dy[j-1], (mesh.dx[i] * mesh.dz[k])
                        )
                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i,j-1,k]
                        AV[n] = as
                    end

                    #North Coefficents
                    if (j != mesh.m1) && (mesh.m1 != 1) && (phi.onoff[i,j+1,k])
                        @inbounds an = _projection_PC_pressure_laplacian_neighbors_(
                            mesh.dy[j], mesh.dy[j+1], (mesh.dx[i] * mesh.dz[k])
                        )
                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i,j+1,k]
                        AV[n] = an
                    end

                    #Bottom Coefficents
                    if (k != 1) && (mesh.n1 != 1) && (phi.onoff[i,j,k-1])
                        @inbounds ab = _projection_PC_pressure_laplacian_neighbors_(
                            mesh.dz[k], mesh.dz[k-1], (mesh.dx[i] * mesh.dy[j])
                        )
                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i,j,k-1]
                        AV[n] = ab
                    end

                    #Top Coefficents
                    if (k != mesh.n1) && (mesh.n1 != 1) && (phi.onoff[i,j,k+1])
                        @inbounds at = _projection_PC_pressure_laplacian_neighbors_(
                            mesh.dz[k], mesh.dz[k+1], (mesh.dx[i] * mesh.dy[j])
                        )
                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i,j,k+1]
                        AV[n] = at
                    end

                    #Center Coefficent
                    n += 1
                    AI[n] = id
                    AJ[n] = id
                    AV[n] = -1.0 * (aw + ae + as + an + ab + at)
                end
            end
        end
    end

    A = sparse(AI[1:n], AJ[1:n], AV[1:n])

    return A
end
