"""

"""
function discretize_SIMPLE_PressureCorrection end

function discretize_SIMPLE_PressureCorrection(
    AU::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,2},
    },
    vel::CSVelocity1D,
    phi::CSPhi1D,
    material::CSMaterial1D,
    mesh::UnionCSMesh1D;
    velocityU::Array{<:AbstractFloat,1} = vel.fValues.uFace,
    T::Type{<:AbstractFloat} = Float64,
    threads::Bool = false,
    sparse::Bool = true,
    interpolation::Signed = 1,
)
    n_equations = maximum_globalIndex(phi)

    if sparse
        A = spzeros(T, n_equations, n_equations)
    elseif !sparse
        A = zeros(T, n_equations, n_equations)
    end

    b = zeros(T, n_equations)

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            if phi.onoff[i]
                id = phi.gIndex[i]

                #Auxiliar variables
                ac = 0.0
                b0 = 0.0
                aw = 0.0
                Dw = 0.0
                b1 = 0.0
                ae = 0.0
                De = 0.0
                b2 = 0.0

                #West Coefficents
                if (i != 1) && (phi.onoff[i-1])
                    rho = density_interpolation(
                        [mesh.dx[i]; mesh.dx[i-1]], [material.ρ[i]; material.ρ[i-1]];
                        interpolation = interpolation
                    )

                    D1 = mesh.vol[i]/AU[phi.gIndex[i],phi.gIndex[i]]
                    D2 = mesh.vol[i-1]/AU[phi.gIndex[i-1],phi.gIndex[i-1]]
                    Dw = general_interpolation(
                        [mesh.dx[i]; mesh.dx[i-1]], [D1; D2];
                        interpolation = interpolation
                    )

                    num = rho * Dw * 1.0
                    den = mesh.x[i] - mesh.x[i-1]
                    aw = num / den

                    b1 = rho * velocityU[i] * 1.0

                    A[id, phi.gIndex[i-1]] = (-1.0) * aw
                    b[id] += b1

                else
                    b1 = material.ρ[i] * velocityU[i] * 1.0

                    b[id] += b1
                end

                #East Coefficents
                if (i != mesh.l1)  && (phi.onoff[i+1])
                    rho = density_interpolation(
                        [mesh.dx[i]; mesh.dx[i+1]], [material.ρ[i]; material.ρ[i+1]];
                        interpolation = interpolation
                    )
                    D1 = mesh.vol[i]/AU[phi.gIndex[i],phi.gIndex[i]]
                    D2 = mesh.vol[i+1]/AU[phi.gIndex[i+1],phi.gIndex[i+1]]
                    De = general_interpolation(
                        [mesh.dx[i]; mesh.dx[i+1]], [D1; D2];
                        interpolation = interpolation
                    )

                    num = rho * De * 1.0
                    den = mesh.x[i+1] - mesh.x[i]
                    ae = num / den

                    b2 = -1.0 * rho * velocityU[i+1] * 1.0

                    A[id, phi.gIndex[i+1]] = (-1.0) * ae
                    b[id] += b2

                else
                    b2 = -1.0 * material.ρ[i] * velocityU[i+1] * 1.0

                    b[id] += b2
                end

                #Center Coefficent
                A[id, id] = ae + aw

            end
        end

    elseif !threads
        for i in 1:mesh.l1
            if phi.onoff[i]
                id = phi.gIndex[i]

                #Auxiliar variables
                ac = 0.0
                b0 = 0.0
                aw = 0.0
                Dw = 0.0
                b1 = 0.0
                ae = 0.0
                De = 0.0
                b2 = 0.0

                #West Coefficents
                if (i != 1) && (phi.onoff[i-1])
                    rho = density_interpolation(
                        [mesh.dx[i]; mesh.dx[i-1]], [material.ρ[i]; material.ρ[i-1]];
                        interpolation = interpolation
                    )

                    D1 = mesh.vol[i]/AU[phi.gIndex[i],phi.gIndex[i]]
                    D2 = mesh.vol[i-1]/AU[phi.gIndex[i-1],phi.gIndex[i-1]]
                    Dw = general_interpolation(
                        [mesh.dx[i]; mesh.dx[i-1]], [D1; D2];
                        interpolation = interpolation
                    )

                    num = rho * Dw * 1.0
                    den = mesh.x[i] - mesh.x[i-1]
                    aw = num / den

                    b1 = rho * velocityU[i] * 1.0

                    A[id, phi.gIndex[i-1]] = (-1.0) * aw
                    b[id] += b1

                else
                    b1 = material.ρ[i] * velocityU[i] * 1.0

                    b[id] += b1
                end

                #East Coefficents
                if (i != mesh.l1)  && (phi.onoff[i+1])
                    rho = density_interpolation(
                        [mesh.dx[i]; mesh.dx[i+1]], [material.ρ[i]; material.ρ[i+1]];
                        interpolation = interpolation
                    )
                    D1 = mesh.vol[i]/AU[phi.gIndex[i],phi.gIndex[i]]
                    D2 = mesh.vol[i+1]/AU[phi.gIndex[i+1],phi.gIndex[i+1]]
                    De = general_interpolation(
                        [mesh.dx[i]; mesh.dx[i+1]], [D1; D2];
                        interpolation = interpolation
                    )

                    num = rho * De * 1.0
                    den = mesh.x[i+1] - mesh.x[i]
                    ae = num / den

                    b2 = -1.0 * rho * velocityU[i+1] * 1.0

                    A[id, phi.gIndex[i+1]] = (-1.0) * ae
                    b[id] += b2

                else
                    b2 = -1.0 * material.ρ[i] * velocityU[i+1] * 1.0

                    b[id] += b2
                end

                #Center Coefficent
                A[id, id] = ae + aw

            end
        end

    end

    if sparse
        A = dropzeros(A)
    end

    return A, b
end

function discretize_SIMPLE_PressureCorrection(
    AU::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,2},
    },
    vel::CSVelocity1D,
    phi::CSPhi1D,
    material::UnionCSConstantMaterial,
    mesh::UnionCSMesh1D;
    velocityU::Array{<:AbstractFloat,1} = vel.fValues.uFace,
    T::Type{<:AbstractFloat} = Float64,
    threads::Bool = false,
    sparse::Bool = true,
    interpolation::Signed = 1,
)
    n_equations = maximum_globalIndex(phi)

    if sparse
        A = spzeros(T, n_equations, n_equations)
    elseif !sparse
        A = zeros(T, n_equations, n_equations)
    end

    b = zeros(T, n_equations)

    if threads
        for i in 1:mesh.l1
            if phi.onoff[i]
                id = phi.gIndex[i]

                #Auxiliar variables
                ac = 0.0
                b0 = 0.0
                aw = 0.0
                Dw = 0.0
                b1 = 0.0
                ae = 0.0
                De = 0.0
                b2 = 0.0

                #West Coefficents
                if (i != 1) && (phi.onoff[i-1])
                    rho = material.ρ

                    D1 = mesh.vol[i]/AU[phi.gIndex[i],phi.gIndex[i]]
                    D2 = mesh.vol[i-1]/AU[phi.gIndex[i-1],phi.gIndex[i-1]]
                    Dw = general_interpolation(
                        [mesh.dx[i]; mesh.dx[i-1]], [D1; D2];
                        interpolation = interpolation
                    )

                    num = rho * Dw * 1.0
                    den = mesh.x[i] - mesh.x[i-1]
                    aw = num / den

                    b1 = rho * velocityU[i] * 1.0

                    A[id, phi.gIndex[i-1]] = (-1.0) * aw
                    b[id] += b1
                else
                    b1 = material.ρ * velocityU[i] * 1.0

                    b[id] += b1
                end

                #East Coefficents
                if (i != mesh.l1)  && (phi.onoff[i+1])
                    rho = material.ρ

                    D1 = mesh.vol[i]/AU[phi.gIndex[i],phi.gIndex[i]]
                    D2 = mesh.vol[i+1]/AU[phi.gIndex[i+1],phi.gIndex[i+1]]
                    De = general_interpolation(
                        [mesh.dx[i]; mesh.dx[i+1]], [D1; D2];
                        interpolation = interpolation
                    )

                    num = rho * De * 1.0
                    den = mesh.x[i+1] - mesh.x[i]
                    ae = num / den

                    b2 = -1.0 * rho * velocityU[i+1] * 1.0

                    A[id, phi.gIndex[i+1]] = (-1.0) * ae
                    b[id] += b2
                else
                    b2 = -1.0 * material.ρ * velocityU[i+1] * 1.0

                    b[id] += b2
                end

                #Center Coefficent
                A[id, id] = ae + aw

            end
        end

    elseif !threads
        Base.Threads.@threads for i in 1:mesh.l1
            if phi.onoff[i]
                id = phi.gIndex[i]

                #Auxiliar variables
                ac = 0.0
                b0 = 0.0
                aw = 0.0
                Dw = 0.0
                b1 = 0.0
                ae = 0.0
                De = 0.0
                b2 = 0.0

                #West Coefficents
                if (i != 1) && (phi.onoff[i-1])
                    rho = material.ρ

                    D1 = mesh.vol[i]/AU[phi.gIndex[i],phi.gIndex[i]]
                    D2 = mesh.vol[i-1]/AU[phi.gIndex[i-1],phi.gIndex[i-1]]
                    Dw = general_interpolation(
                        [mesh.dx[i]; mesh.dx[i-1]], [D1; D2];
                        interpolation = interpolation
                    )

                    num = rho * Dw * 1.0
                    den = mesh.x[i] - mesh.x[i-1]
                    aw = num / den

                    b1 = rho * velocityU[i] * 1.0

                    A[id, phi.gIndex[i-1]] = (-1.0) * aw
                    b[id] += b1
                else
                    b1 = material.ρ * velocityU[i] * 1.0

                    b[id] += b1
                end

                #East Coefficents
                if (i != mesh.l1)  && (phi.onoff[i+1])
                    rho = material.ρ

                    D1 = mesh.vol[i]/AU[phi.gIndex[i],phi.gIndex[i]]
                    D2 = mesh.vol[i+1]/AU[phi.gIndex[i+1],phi.gIndex[i+1]]
                    De = general_interpolation(
                        [mesh.dx[i]; mesh.dx[i+1]], [D1; D2];
                        interpolation = interpolation
                    )

                    num = rho * De * 1.0
                    den = mesh.x[i+1] - mesh.x[i]
                    ae = num / den

                    b2 = -1.0 * rho * velocityU[i+1] * 1.0

                    A[id, phi.gIndex[i+1]] = (-1.0) * ae
                    b[id] += b2
                else
                    b2 = -1.0 * material.ρ * velocityU[i+1] * 1.0

                    b[id] += b2
                end

                #Center Coefficent
                A[id, id] = ae + aw

            end
        end

    end

    if sparse
        A = dropzeros(A)
    end

    return A, b
end

function discretize_SIMPLE_PressureCorrection(
    AU::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,2},
    },
    AV::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,2},
    },
    vel::CSVelocity2D,
    phi::CSPhi2D,
    material::CSMaterial2D,
    mesh::UnionCSMesh2D;
    velocityU::Array{<:AbstractFloat,2} = vel.fValues.uFace,
    velocityV::Array{<:AbstractFloat,2} = vel.fValues.vFace,
    T::Type{<:AbstractFloat} = Float64,
    threads::Bool = false,
    sparse::Bool = true,
    scheme::Signed = 2,
    interpolation::Signed = 1,
)
    n_equations = maximum_globalIndex(phi)

    if sparse
        A = spzeros(T, n_equations, n_equations)
    elseif !sparse
        A = zeros(T, n_equations, n_equations)
    end

    b = zeros(T, n_equations)

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            for j in 1:mesh.m1
                if phi.onoff[i,j]
                    id = phi.gIndex[i,j]

                    #Auxiliar variables
                    ac = 0.0
                    b0 = 0.0
                    aw = 0.0
                    Dw = 0.0
                    b1 = 0.0
                    ae = 0.0
                    De = 0.0
                    b2 = 0.0
                    as = 0.0
                    Ds = 0.0
                    b3 = 0.0
                    an = 0.0
                    Dn = 0.0
                    b4 = 0.0

                    #West Coefficents
                    if (i != 1) && (phi.onoff[i-1,j])
                        rho = density_interpolation(
                            [mesh.dx[i]; mesh.dx[i-1]], [material.ρ[i,j]; material.ρ[i-1,j]];
                            interpolation = interpolation
                        )

                        D1 = mesh.vol[i,j]/AU[phi.gIndex[i,j],phi.gIndex[i,j]]
                        D2 = mesh.vol[i-1,j]/AU[phi.gIndex[i-1,j],phi.gIndex[i-1,j]]
                        Dw = general_interpolation(
                            [mesh.dx[i]; mesh.dx[i-1]], [D1; D2];
                            interpolation = interpolation
                        )

                        num = rho * Dw * mesh.dy[j]
                        den = mesh.x[i] - mesh.x[i-1]
                        aw = num / den

                        b1 = rho * velocityU[i,j] * mesh.dy[j]

                        A[id, phi.gIndex[i-1,j]] = (-1.0) * aw
                        b[id] += b1
                    else
                        b1 = material.ρ[i,j] * velocityU[i,j] * mesh.dy[j]

                        b[id] += b1
                    end

                    #East Coefficents
                    if (i != mesh.l1)  && (phi.onoff[i+1,j])
                        rho = density_interpolation(
                            [mesh.dx[i]; mesh.dx[i+1]], [material.ρ[i,j]; material.ρ[i+1,j]];
                            interpolation = interpolation
                        )
                        D1 = mesh.vol[i,j]/AU[phi.gIndex[i,j],phi.gIndex[i,j]]
                        D2 = mesh.vol[i+1,j]/AU[phi.gIndex[i+1,j],phi.gIndex[i+1,j]]
                        De = general_interpolation(
                            [mesh.dx[i]; mesh.dx[i+1]], [D1; D2];
                            interpolation = interpolation
                        )

                        num = rho * De * mesh.dy[j]
                        den = mesh.x[i+1] - mesh.x[i]
                        ae = num / den

                        b2 = -1.0 * rho * velocityU[i+1,j] * mesh.dy[j]

                        A[id, phi.gIndex[i+1,j]] = (-1.0) * ae
                        b[id] += b2
                    else
                        b2 = -1.0 * material.ρ[i,j] * velocityU[i+1,j] * mesh.dy[j]

                        b[id] += b2
                    end

                    #South Coefficents
                    if (j != 1) && (phi.onoff[i,j-1])
                        rho = density_interpolation(
                            [mesh.dy[j]; mesh.dy[j-1]], [material.ρ[i,j]; material.ρ[i,j-1]];
                            interpolation = interpolation
                        )

                        D1 = mesh.vol[i,j]/AV[phi.gIndex[i,j],phi.gIndex[i,j]]
                        D2 = mesh.vol[i,j-1]/AV[phi.gIndex[i,j-1],phi.gIndex[i,j-1]]
                        Ds = general_interpolation(
                            [mesh.dy[j]; mesh.dy[j-1]], [D1; D2];
                            interpolation = interpolation
                        )

                        num = rho * Ds * mesh.dx[i]
                        den = mesh.y[j] - mesh.y[j-1]
                        as = num / den

                        b3 = rho * velocityV[i,j] * mesh.dx[i]

                        A[id, phi.gIndex[i,j-1]] = (-1.0) * as
                        b[id] += b3
                    else
                        b3 = material.ρ[i,j] * velocityV[i,j] * mesh.dx[i]

                        b[id] += b3
                    end

                    #North Coefficents
                    if (j != mesh.m1) && (phi.onoff[i,j+1])
                        rho = density_interpolation(
                            [mesh.dy[j]; mesh.dy[j+1]], [material.ρ[i,j]; material.ρ[i,j+1]];
                            interpolation = interpolation
                        )
                        D1 = mesh.vol[i,j]/AV[phi.gIndex[i,j],phi.gIndex[i,j]]
                        D2 = mesh.vol[i,j+1]/AV[phi.gIndex[i,j+1],phi.gIndex[i,j+1]]
                        Dn = general_interpolation(
                            [mesh.dy[j]; mesh.dy[j+1]], [D1; D2];
                            interpolation = interpolation
                        )

                        num = rho * Dn * mesh.dx[i]
                        den = mesh.y[j+1] - mesh.y[j]
                        an = num / den

                        b4 = -1.0 * rho * velocityV[i,j+1] * mesh.dx[i]

                        A[id, phi.gIndex[i,j+1]] = (-1.0) * an
                        b[id] += b4
                    else
                        b4 = -1.0 * material.ρ[i,j] * velocityV[i,j+1] * mesh.dx[i]

                        b[id] += b4
                    end

                    #Center Coefficent
                    A[id, id] = ae + aw + as + an

                end
            end
        end

    elseif !threads
        for i in 1:mesh.l1
            for j in 1:mesh.m1
                if phi.onoff[i,j]
                    id = phi.gIndex[i,j]

                    #Auxiliar variables
                    ac = 0.0
                    b0 = 0.0
                    aw = 0.0
                    Dw = 0.0
                    b1 = 0.0
                    ae = 0.0
                    De = 0.0
                    b2 = 0.0
                    as = 0.0
                    Ds = 0.0
                    b3 = 0.0
                    an = 0.0
                    Dn = 0.0
                    b4 = 0.0

                    #West Coefficents
                    if (i != 1) && (phi.onoff[i-1,j])
                        rho = density_interpolation(
                            [mesh.dx[i]; mesh.dx[i-1]], [material.ρ[i,j]; material.ρ[i-1,j]];
                            interpolation = interpolation
                        )

                        D1 = mesh.vol[i,j]/AU[phi.gIndex[i,j],phi.gIndex[i,j]]
                        D2 = mesh.vol[i-1,j]/AU[phi.gIndex[i-1,j],phi.gIndex[i-1,j]]
                        Dw = general_interpolation(
                            [mesh.dx[i]; mesh.dx[i-1]], [D1; D2];
                            interpolation = interpolation
                        )

                        num = rho * Dw * mesh.dy[j]
                        den = mesh.x[i] - mesh.x[i-1]
                        aw = num / den

                        b1 = rho * velocityU[i,j] * mesh.dy[j]

                        A[id, phi.gIndex[i-1,j]] = (-1.0) * aw
                        b[id] += b1
                    else
                        b1 = material.ρ[i,j] * velocityU[i,j] * mesh.dy[j]

                        b[id] += b1
                    end

                    #East Coefficents
                    if (i != mesh.l1)  && (phi.onoff[i+1,j])
                        rho = density_interpolation(
                            [mesh.dx[i]; mesh.dx[i+1]], [material.ρ[i,j]; material.ρ[i+1,j]];
                            interpolation = interpolation
                        )
                        D1 = mesh.vol[i,j]/AU[phi.gIndex[i,j],phi.gIndex[i,j]]
                        D2 = mesh.vol[i+1,j]/AU[phi.gIndex[i+1,j],phi.gIndex[i+1,j]]
                        De = general_interpolation(
                            [mesh.dx[i]; mesh.dx[i+1]], [D1; D2];
                            interpolation = interpolation
                        )

                        num = rho * De * mesh.dy[j]
                        den = mesh.x[i+1] - mesh.x[i]
                        ae = num / den

                        b2 = -1.0 * rho * velocityU[i+1,j] * mesh.dy[j]

                        A[id, phi.gIndex[i+1,j]] = (-1.0) * ae
                        b[id] += b2
                    else
                        b2 = -1.0 * material.ρ[i,j] * velocityU[i+1,j] * mesh.dy[j]

                        b[id] += b2
                    end

                    #South Coefficents
                    if (j != 1) && (phi.onoff[i,j-1])
                        rho = density_interpolation(
                            [mesh.dy[j]; mesh.dy[j-1]], [material.ρ[i,j]; material.ρ[i,j-1]];
                            interpolation = interpolation
                        )

                        D1 = mesh.vol[i,j]/AV[phi.gIndex[i,j],phi.gIndex[i,j]]
                        D2 = mesh.vol[i,j-1]/AV[phi.gIndex[i,j-1],phi.gIndex[i,j-1]]
                        Ds = general_interpolation(
                            [mesh.dy[j]; mesh.dy[j-1]], [D1; D2];
                            interpolation = interpolation
                        )

                        num = rho * Ds * mesh.dx[i]
                        den = mesh.y[j] - mesh.y[j-1]
                        as = num / den

                        b3 = rho * velocityV[i,j] * mesh.dx[i]

                        A[id, phi.gIndex[i,j-1]] = (-1.0) * as
                        b[id] += b3
                    else
                        b3 = material.ρ[i,j] * velocityV[i,j] * mesh.dx[i]

                        b[id] += b3
                    end

                    #North Coefficents
                    if (j != mesh.m1) && (phi.onoff[i,j+1])
                        rho = density_interpolation(
                            [mesh.dy[j]; mesh.dy[j+1]], [material.ρ[i,j]; material.ρ[i,j+1]];
                            interpolation = interpolation
                        )
                        D1 = mesh.vol[i,j]/AV[phi.gIndex[i,j],phi.gIndex[i,j]]
                        D2 = mesh.vol[i,j+1]/AV[phi.gIndex[i,j+1],phi.gIndex[i,j+1]]
                        Dn = general_interpolation(
                            [mesh.dy[j]; mesh.dy[j+1]], [D1; D2];
                            interpolation = interpolation
                        )

                        num = rho * Dn * mesh.dx[i]
                        den = mesh.y[j+1] - mesh.y[j]
                        an = num / den

                        b4 = -1.0 * rho * velocityV[i,j+1] * mesh.dx[i]

                        A[id, phi.gIndex[i,j+1]] = (-1.0) * an
                        b[id] += b4
                    else
                        b4 = -1.0 * material.ρ[i,j] * velocityV[i,j+1] * mesh.dx[i]

                        b[id] += b4
                    end

                    #Center Coefficent
                    A[id, id] = ae + aw + as + an

                end
            end
        end

    end

    if sparse
        A = dropzeros(A)
    end

    return A, b
end

function discretize_SIMPLE_PressureCorrection(
    AU::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,2},
    },
    AV::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,2},
    },
    vel::CSVelocity2D,
    phi::CSPhi2D,
    material::UnionCSConstantMaterial,
    mesh::UnionCSMesh2D;
    velocityU::Array{<:AbstractFloat,2} = vel.fValues.uFace,
    velocityV::Array{<:AbstractFloat,2} = vel.fValues.vFace,
    T::Type{<:AbstractFloat} = Float64,
    threads::Bool = false,
    sparse::Bool = true,
    scheme::Signed = 2,
    interpolation::Signed = 1,
)
    n_equations = maximum_globalIndex(phi)

    if sparse
        A = spzeros(T, n_equations, n_equations)
    elseif !sparse
        A = zeros(T, n_equations, n_equations)
    end

    b = zeros(T, n_equations)

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            for j in 1:mesh.m1
                if phi.onoff[i,j]
                    id = phi.gIndex[i,j]

                    #Auxiliar variables
                    ac = 0.0
                    b0 = 0.0
                    aw = 0.0
                    Dw = 0.0
                    b1 = 0.0
                    ae = 0.0
                    De = 0.0
                    b2 = 0.0
                    as = 0.0
                    Ds = 0.0
                    b3 = 0.0
                    an = 0.0
                    Dn = 0.0
                    b4 = 0.0

                    #West Coefficents
                    if (i != 1) && (phi.onoff[i-1,j])
                        rho = material.ρ

                        D1 = mesh.vol[i,j]/AU[phi.gIndex[i,j],phi.gIndex[i,j]]
                        D2 = mesh.vol[i-1,j]/AU[phi.gIndex[i-1,j],phi.gIndex[i-1,j]]
                        Dw = general_interpolation(
                            [mesh.dx[i]; mesh.dx[i-1]], [D1; D2];
                            interpolation = interpolation
                        )

                        num = rho * Dw * mesh.dy[j]
                        den = mesh.x[i] - mesh.x[i-1]
                        aw = num / den

                        b1 = rho * velocityU[i,j] * mesh.dy[j]

                        A[id, phi.gIndex[i-1,j]] = (-1.0) * aw
                        b[id] += b1
                    else
                        b1 = material.ρ * velocityU[i,j] * mesh.dy[j]

                        b[id] += b1
                    end

                    #East Coefficents
                    if (i != mesh.l1)  && (phi.onoff[i+1,j])
                        rho = material.ρ

                        D1 = mesh.vol[i,j]/AU[phi.gIndex[i,j],phi.gIndex[i,j]]
                        D2 = mesh.vol[i+1,j]/AU[phi.gIndex[i+1,j],phi.gIndex[i+1,j]]
                        De = general_interpolation(
                            [mesh.dx[i]; mesh.dx[i+1]], [D1; D2];
                            interpolation = interpolation
                        )

                        num = rho * De * mesh.dy[j]
                        den = mesh.x[i+1] - mesh.x[i]
                        ae = num / den

                        b2 = -1.0 * rho * velocityU[i+1,j] * mesh.dy[j]

                        A[id, phi.gIndex[i+1,j]] = (-1.0) * ae
                        b[id] += b2
                    else
                        b2 = -1.0 * material.ρ * velocityU[i+1,j] * mesh.dy[j]

                        b[id] += b2
                    end

                    #South Coefficents
                    if (j != 1) && (phi.onoff[i,j-1])
                        rho = material.ρ

                        D1 = mesh.vol[i,j]/AV[phi.gIndex[i,j],phi.gIndex[i,j]]
                        D2 = mesh.vol[i,j-1]/AV[phi.gIndex[i,j-1],phi.gIndex[i,j-1]]
                        Ds = general_interpolation(
                            [mesh.dy[j]; mesh.dy[j-1]], [D1; D2];
                            interpolation = interpolation
                        )

                        num = rho * Ds * mesh.dx[i]
                        den = mesh.y[j] - mesh.y[j-1]
                        as = num / den

                        b3 = rho * velocityV[i,j] * mesh.dx[i]

                        A[id, phi.gIndex[i,j-1]] = (-1.0) * as
                        b[id] += b3
                    else
                        b3 = material.ρ * velocityV[i,j] * mesh.dx[i]

                        b[id] += b3
                    end

                    #North Coefficents
                    if (j != mesh.m1) && (phi.onoff[i,j+1])
                        rho = material.ρ

                        D1 = mesh.vol[i,j]/AV[phi.gIndex[i,j],phi.gIndex[i,j]]
                        D2 = mesh.vol[i,j+1]/AV[phi.gIndex[i,j+1],phi.gIndex[i,j+1]]
                        Dn = general_interpolation(
                            [mesh.dy[j]; mesh.dy[j+1]], [D1; D2];
                            interpolation = interpolation
                        )

                        num = rho * Dn * mesh.dx[i]
                        den = mesh.y[j+1] - mesh.y[j]
                        an = num / den

                        b4 = -1.0 * rho * velocityV[i,j+1] * mesh.dx[i]

                        A[id, phi.gIndex[i,j+1]] = (-1.0) * an
                        b[id] += b4
                    else
                        b4 = -1.0 * material.ρ * velocityV[i,j+1] * mesh.dx[i]

                        b[id] += b4
                    end

                    #Center Coefficent
                    A[id, id] = ae + aw + as + an

                end
            end
        end

    elseif !threads
        for i in 1:mesh.l1
            for j in 1:mesh.m1
                if phi.onoff[i,j]
                    id = phi.gIndex[i,j]

                    #Auxiliar variables
                    ac = 0.0
                    b0 = 0.0
                    aw = 0.0
                    Dw = 0.0
                    b1 = 0.0
                    ae = 0.0
                    De = 0.0
                    b2 = 0.0
                    as = 0.0
                    Ds = 0.0
                    b3 = 0.0
                    an = 0.0
                    Dn = 0.0
                    b4 = 0.0

                    #West Coefficents
                    if (i != 1) && (phi.onoff[i-1,j])
                        rho = material.ρ

                        D1 = mesh.vol[i,j]/AU[phi.gIndex[i,j],phi.gIndex[i,j]]
                        D2 = mesh.vol[i-1,j]/AU[phi.gIndex[i-1,j],phi.gIndex[i-1,j]]
                        Dw = general_interpolation(
                            [mesh.dx[i]; mesh.dx[i-1]], [D1; D2];
                            interpolation = interpolation
                        )

                        num = rho * Dw * mesh.dy[j]
                        den = mesh.x[i] - mesh.x[i-1]
                        aw = num / den

                        b1 = rho * velocityU[i,j] * mesh.dy[j]

                        A[id, phi.gIndex[i-1,j]] = (-1.0) * aw
                        b[id] += b1
                    else
                        b1 = material.ρ * velocityU[i,j] * mesh.dy[j]

                        b[id] += b1
                    end

                    #East Coefficents
                    if (i != mesh.l1)  && (phi.onoff[i+1,j])
                        rho = material.ρ

                        D1 = mesh.vol[i,j]/AU[phi.gIndex[i,j],phi.gIndex[i,j]]
                        D2 = mesh.vol[i+1,j]/AU[phi.gIndex[i+1,j],phi.gIndex[i+1,j]]
                        De = general_interpolation(
                            [mesh.dx[i]; mesh.dx[i+1]], [D1; D2];
                            interpolation = interpolation
                        )

                        num = rho * De * mesh.dy[j]
                        den = mesh.x[i+1] - mesh.x[i]
                        ae = num / den

                        b2 = -1.0 * rho * velocityU[i+1,j] * mesh.dy[j]

                        A[id, phi.gIndex[i+1,j]] = (-1.0) * ae
                        b[id] += b2
                    else
                        b2 = -1.0 * material.ρ * velocityU[i+1,j] * mesh.dy[j]

                        b[id] += b2
                    end

                    #South Coefficents
                    if (j != 1) && (phi.onoff[i,j-1])
                        rho = material.ρ

                        D1 = mesh.vol[i,j]/AV[phi.gIndex[i,j],phi.gIndex[i,j]]
                        D2 = mesh.vol[i,j-1]/AV[phi.gIndex[i,j-1],phi.gIndex[i,j-1]]
                        Ds = general_interpolation(
                            [mesh.dy[j]; mesh.dy[j-1]], [D1; D2];
                            interpolation = interpolation
                        )

                        num = rho * Ds * mesh.dx[i]
                        den = mesh.y[j] - mesh.y[j-1]
                        as = num / den

                        b3 = rho * velocityV[i,j] * mesh.dx[i]

                        A[id, phi.gIndex[i,j-1]] = (-1.0) * as
                        b[id] += b3
                    else
                        b3 = material.ρ * velocityV[i,j] * mesh.dx[i]

                        b[id] += b3
                    end

                    #North Coefficents
                    if (j != mesh.m1) && (phi.onoff[i,j+1])
                        rho = material.ρ

                        D1 = mesh.vol[i,j]/AV[phi.gIndex[i,j],phi.gIndex[i,j]]
                        D2 = mesh.vol[i,j+1]/AV[phi.gIndex[i,j+1],phi.gIndex[i,j+1]]
                        Dn = general_interpolation(
                            [mesh.dy[j]; mesh.dy[j+1]], [D1; D2];
                            interpolation = interpolation
                        )

                        num = rho * Dn * mesh.dx[i]
                        den = mesh.y[j+1] - mesh.y[j]
                        an = num / den

                        b4 = -1.0 * rho * velocityV[i,j+1] * mesh.dx[i]

                        A[id, phi.gIndex[i,j+1]] = (-1.0) * an
                        b[id] += b4
                    else
                        b4 = -1.0 * material.ρ * velocityV[i,j+1] * mesh.dx[i]

                        b[id] += b4
                    end

                    #Center Coefficent
                    A[id, id] = ae + aw + as + an

                end
            end
        end

    end

    if sparse
        A = dropzeros(A)
    end

    return A, b
end


function discretize_SIMPLE_PressureCorrection(
    AU::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,2},
    },
    AV::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,2},
    },
    AW::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,2},
    },
    vel::CSVelocity3D,
    phi::CSPhi3D,
    material::CSMaterial3D,
    mesh::UnionCSMesh3D;
    velocityU::Array{<:AbstractFloat,3} = vel.fValues.uFace,
    velocityV::Array{<:AbstractFloat,3} = vel.fValues.vFace,
    velocityW::Array{<:AbstractFloat,3} = vel.fValues.wFace,
    T::Type{<:AbstractFloat} = Float64,
    threads::Bool = false,
    sparse::Bool = true,
    scheme::Signed = 2,
    interpolation::Signed = 1,
)
    n_equations = maximum_globalIndex(phi)

    if sparse
        A = spzeros(T, n_equations, n_equations)
    elseif !sparse
        A = zeros(T, n_equations, n_equations)
    end

    b = zeros(T, n_equations)

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            for j in 1:mesh.m1
                for k in 1:mesh.n1
                    if phi.onoff[i,j,k]
                        id = phi.gIndex[i,j,k]

                        #Auxiliar variables
                        ac = 0.0
                        b0 = 0.0
                        aw = 0.0
                        Dw = 0.0
                        b1 = 0.0
                        ae = 0.0
                        De = 0.0
                        b2 = 0.0
                        as = 0.0
                        Ds = 0.0
                        b3 = 0.0
                        an = 0.0
                        Dn = 0.0
                        b4 = 0.0
                        ab = 0.0
                        Db = 0.0
                        b5 = 0.0
                        at = 0.0
                        Dt = 0.0
                        b6 = 0.0

                        #West Coefficents
                        if (i != 1) && (phi.onoff[i-1,j,k])
                            rho = density_interpolation(
                                [mesh.dx[i]; mesh.dx[i-1]], [material.ρ[i,j,k]; material.ρ[i-1,j,k]];
                                interpolation = interpolation
                            )

                            D1 = mesh.vol[i,j,k]/AU[phi.gIndex[i,j,k],phi.gIndex[i,j,k]]
                            D2 = mesh.vol[i-1,j,k]/AU[phi.gIndex[i-1,j,k],phi.gIndex[i-1,j,k]]
                            Dw = general_interpolation(
                                [mesh.dx[i]; mesh.dx[i-1]], [D1; D2];
                                interpolation = interpolation
                            )

                            num = rho * Dw * (mesh.dy[j] * mesh.dz[k])
                            den = mesh.x[i] - mesh.x[i-1]
                            aw = num / den

                            b1 = rho * velocityU[i,j,k] * (mesh.dy[j] * mesh.dz[k])

                            A[id, phi.gIndex[i-1,j,k]] = (-1.0) * aw
                            b[id] += b1
                        else
                            b1 = material.ρ[i,j,k] * velocityU[i,j,k] * (mesh.dy[j] * mesh.dz[k])

                            b[id] += b1
                        end

                        #East Coefficents
                        if (i != mesh.l1)  && (phi.onoff[i+1,j,k])
                            rho = density_interpolation(
                                [mesh.dx[i]; mesh.dx[i+1]], [material.ρ[i,j,k]; material.ρ[i+1,j,k]];
                                interpolation = interpolation
                            )
                            D1 = mesh.vol[i,j,k]/AU[phi.gIndex[i,j,k],phi.gIndex[i,j,k]]
                            D2 = mesh.vol[i+1,j,k]/AU[phi.gIndex[i+1,j,k],phi.gIndex[i+1,j,k]]
                            De = general_interpolation(
                                [mesh.dx[i]; mesh.dx[i+1]], [D1; D2];
                                interpolation = interpolation
                            )

                            num = rho * De * (mesh.dy[j] * mesh.dz[k])
                            den = mesh.x[i+1] - mesh.x[i]
                            ae = num / den

                            b2 = -1.0 * rho * velocityU[i+1,j,k] * (mesh.dy[j] * mesh.dz[k])

                            A[id, phi.gIndex[i+1,j,k]] = (-1.0) * ae
                            b[id] += b2
                        else
                            b2 = -1.0 * material.ρ[i,j,k] * velocityU[i+1,j,k] * (mesh.dy[j] * mesh.dz[k])

                            b[id] += b2
                        end

                        #South Coefficents
                        if (j != 1) && (phi.onoff[i,j-1,k])
                            rho = density_interpolation(
                                [mesh.dy[j]; mesh.dy[j-1]], [material.ρ[i,j,k]; material.ρ[i,j-1,k]];
                                interpolation = interpolation
                            )

                            D1 = mesh.vol[i,j,k]/AV[phi.gIndex[i,j,k],phi.gIndex[i,j,k]]
                            D2 = mesh.vol[i,j-1,k]/AV[phi.gIndex[i,j-1,k],phi.gIndex[i,j-1,k]]
                            Ds = general_interpolation(
                                [mesh.dy[j]; mesh.dy[j-1]], [D1; D2];
                                interpolation = interpolation
                            )

                            num = rho * Ds * (mesh.dx[i] * mesh.dz[k])
                            den = mesh.y[j] - mesh.y[j-1]
                            as = num / den

                            b3 = rho * velocityV[i,j,k] * (mesh.dx[i] * mesh.dz[k])

                            A[id, phi.gIndex[i,j-1,k]] = (-1.0) * as
                            b[id] += b3
                        else
                            b3 = material.ρ[i,j,k] * velocityV[i,j,k] * (mesh.dx[i] * mesh.dz[k])

                            b[id] += b3
                        end

                        #North Coefficents
                        if (j != mesh.l1) && (phi.onoff[i,j+1,k])
                            rho = density_interpolation(
                                [mesh.dy[j]; mesh.dy[j+1]], [material.ρ[i,j,k]; material.ρ[i,j+1,k]];
                                interpolation = interpolation
                            )
                            D1 = mesh.vol[i,j,k]/AV[phi.gIndex[i,j,k],phi.gIndex[i,j,k]]
                            D2 = mesh.vol[i,j+1,k]/AV[phi.gIndex[i,j+1,k],phi.gIndex[i,j+1,k]]
                            Dn = general_interpolation(
                                [mesh.dy[j]; mesh.dy[j+1]], [D1; D2];
                                interpolation = interpolation
                            )

                            num = rho * Dn * (mesh.dx[i] * mesh.dz[k])
                            den = mesh.y[j+1] - mesh.y[j]
                            an = num / den

                            b4 = -1.0 * rho * velocityV[i,j+1,k] * (mesh.dx[i] * mesh.dz[k])

                            A[id, phi.gIndex[i,j+1,k]] = (-1.0) * an
                            b[id] += b4
                        else
                            b4 = -1.0 * material.ρ[i,j,k] * velocityV[i,j+1,k] * (mesh.dx[i] * mesh.dz[k])

                            b[id] += b4
                        end

                        #Bottom Coefficents
                        if (k != 1) && (phi.onoff[i,j,k-1])
                            rho = density_interpolation(
                                [mesh.dz[k]; mesh.dz[k-1]], [material.ρ[i,j,k]; material.ρ[i,j,k-1]];
                                interpolation = interpolation
                            )

                            D1 = mesh.vol[i,j,k]/AW[phi.gIndex[i,j,k],phi.gIndex[i,j,k]]
                            D2 = mesh.vol[i,j,k-1]/AW[phi.gIndex[i,j,k-1],phi.gIndex[i,j,k-1]]
                            Db = general_interpolation(
                                [mesh.dz[k]; mesh.dz[k-1]], [D1; D2];
                                interpolation = interpolation
                            )

                            num = rho * Db * (mesh.dx[i] * mesh.dy[j])
                            den = mesh.z[k] - mesh.z[k-1]
                            ab = num / den

                            b5 = rho * velocityW[i,j,k] * (mesh.dx[i] * mesh.dy[j])

                            A[id, phi.gIndex[i,j,k-1]] = (-1.0) * ab
                            b[id] += b5
                        else
                            b5 = material.ρ[i,j,k] * velocityW[i,j,k] * (mesh.dx[i] * mesh.dy[j])

                            b[id] += b5
                        end

                        #Top Coefficents
                        if (k != mesh.n1) && (phi.onoff[i,j,k+1])
                            rho = density_interpolation(
                                [mesh.dz[k]; mesh.dz[k+1]], [material.ρ[i,j,k]; material.ρ[i,j,k+1]];
                                interpolation = interpolation
                            )
                            D1 = mesh.vol[i,j,k]/AW[phi.gIndex[i,j,k],phi.gIndex[i,j,k]]
                            D2 = mesh.vol[i,j,k+1]/AW[phi.gIndex[i,j,k+1],phi.gIndex[i,j,k+1]]
                            Dt = general_interpolation(
                                [mesh.dz[k]; mesh.dz[k+1]], [D1; D2];
                                interpolation = interpolation
                            )

                            num = rho * Dt * (mesh.dx[i] * mesh.dy[j])
                            den = mesh.z[k+1] - mesh.z[k]
                            at = num / den

                            b6 = -1.0 * rho * velocityW[i,j,k+1] * (mesh.dx[i] * mesh.dy[j])

                            A[id, phi.gIndex[i,j,k+1]] = (-1.0) * at
                            b[id] += b6
                        else
                            b6 = -1.0 * material.ρ[i,j,k] * velocityW[i,j,k+1] * (mesh.dx[i] * mesh.dy[j])

                            b[id] += b6
                        end

                        #Center Coefficent
                        A[id, id] = aw + ae + as + an + ab + at

                    end
                end
            end
        end

    elseif !threads
        for i in 1:mesh.l1
            for j in 1:mesh.m1
                for k in 1:mesh.n1
                    if phi.onoff[i,j,k]
                        id = phi.gIndex[i,j,k]

                        #Auxiliar variables
                        ac = 0.0
                        b0 = 0.0
                        aw = 0.0
                        Dw = 0.0
                        b1 = 0.0
                        ae = 0.0
                        De = 0.0
                        b2 = 0.0
                        as = 0.0
                        Ds = 0.0
                        b3 = 0.0
                        an = 0.0
                        Dn = 0.0
                        b4 = 0.0
                        ab = 0.0
                        Db = 0.0
                        b5 = 0.0
                        at = 0.0
                        Dt = 0.0
                        b6 = 0.0

                        #West Coefficents
                        if (i != 1) && (phi.onoff[i-1,j,k])
                            rho = density_interpolation(
                                [mesh.dx[i]; mesh.dx[i-1]], [material.ρ[i,j,k]; material.ρ[i-1,j,k]];
                                interpolation = interpolation
                            )

                            D1 = mesh.vol[i,j,k]/AU[phi.gIndex[i,j,k],phi.gIndex[i,j,k]]
                            D2 = mesh.vol[i-1,j,k]/AU[phi.gIndex[i-1,j,k],phi.gIndex[i-1,j,k]]
                            Dw = general_interpolation(
                                [mesh.dx[i]; mesh.dx[i-1]], [D1; D2];
                                interpolation = interpolation
                            )

                            num = rho * Dw * (mesh.dy[j] * mesh.dz[k])
                            den = mesh.x[i] - mesh.x[i-1]
                            aw = num / den

                            b1 = rho * velocityU[i,j,k] * (mesh.dy[j] * mesh.dz[k])

                            A[id, phi.gIndex[i-1,j,k]] = (-1.0) * aw
                            b[id] += b1
                        else
                            b1 = material.ρ[i,j,k] * velocityU[i,j,k] * (mesh.dy[j] * mesh.dz[k])

                            b[id] += b1
                        end

                        #East Coefficents
                        if (i != mesh.l1)  && (phi.onoff[i+1,j,k])
                            rho = density_interpolation(
                                [mesh.dx[i]; mesh.dx[i+1]], [material.ρ[i,j,k]; material.ρ[i+1,j,k]];
                                interpolation = interpolation
                            )
                            D1 = mesh.vol[i,j,k]/AU[phi.gIndex[i,j,k],phi.gIndex[i,j,k]]
                            D2 = mesh.vol[i+1,j,k]/AU[phi.gIndex[i+1,j,k],phi.gIndex[i+1,j,k]]
                            De = general_interpolation(
                                [mesh.dx[i]; mesh.dx[i+1]], [D1; D2];
                                interpolation = interpolation
                            )

                            num = rho * De * (mesh.dy[j] * mesh.dz[k])
                            den = mesh.x[i+1] - mesh.x[i]
                            ae = num / den

                            b2 = -1.0 * rho * velocityU[i+1,j,k] * (mesh.dy[j] * mesh.dz[k])

                            A[id, phi.gIndex[i+1,j,k]] = (-1.0) * ae
                            b[id] += b2
                        else
                            b2 = -1.0 * material.ρ[i,j,k] * velocityU[i+1,j,k] * (mesh.dy[j] * mesh.dz[k])

                            b[id] += b2
                        end

                        #South Coefficents
                        if (j != 1) && (phi.onoff[i,j-1,k])
                            rho = density_interpolation(
                                [mesh.dy[j]; mesh.dy[j-1]], [material.ρ[i,j,k]; material.ρ[i,j-1,k]];
                                interpolation = interpolation
                            )

                            D1 = mesh.vol[i,j,k]/AV[phi.gIndex[i,j,k],phi.gIndex[i,j,k]]
                            D2 = mesh.vol[i,j-1,k]/AV[phi.gIndex[i,j-1,k],phi.gIndex[i,j-1,k]]
                            Ds = general_interpolation(
                                [mesh.dy[j]; mesh.dy[j-1]], [D1; D2];
                                interpolation = interpolation
                            )

                            num = rho * Ds * (mesh.dx[i] * mesh.dz[k])
                            den = mesh.y[j] - mesh.y[j-1]
                            as = num / den

                            b3 = rho * velocityV[i,j,k] * (mesh.dx[i] * mesh.dz[k])

                            A[id, phi.gIndex[i,j-1,k]] = (-1.0) * as
                            b[id] += b3
                        else
                            b3 = material.ρ[i,j,k] * velocityV[i,j,k] * (mesh.dx[i] * mesh.dz[k])

                            b[id] += b3
                        end

                        #North Coefficents
                        if (j != mesh.l1) && (phi.onoff[i,j+1,k])
                            rho = density_interpolation(
                                [mesh.dy[j]; mesh.dy[j+1]], [material.ρ[i,j,k]; material.ρ[i,j+1,k]];
                                interpolation = interpolation
                            )
                            D1 = mesh.vol[i,j,k]/AV[phi.gIndex[i,j,k],phi.gIndex[i,j,k]]
                            D2 = mesh.vol[i,j+1,k]/AV[phi.gIndex[i,j+1,k],phi.gIndex[i,j+1,k]]
                            Dn = general_interpolation(
                                [mesh.dy[j]; mesh.dy[j+1]], [D1; D2];
                                interpolation = interpolation
                            )

                            num = rho * Dn * (mesh.dx[i] * mesh.dz[k])
                            den = mesh.y[j+1] - mesh.y[j]
                            an = num / den

                            b4 = -1.0 * rho * velocityV[i,j+1,k] * (mesh.dx[i] * mesh.dz[k])

                            A[id, phi.gIndex[i,j+1,k]] = (-1.0) * an
                            b[id] += b4
                        else
                            b4 = -1.0 * material.ρ[i,j,k] * velocityV[i,j+1,k] * (mesh.dx[i] * mesh.dz[k])

                            b[id] += b4
                        end

                        #Bottom Coefficents
                        if (k != 1) && (phi.onoff[i,j,k-1])
                            rho = density_interpolation(
                                [mesh.dz[k]; mesh.dz[k-1]], [material.ρ[i,j,k]; material.ρ[i,j,k-1]];
                                interpolation = interpolation
                            )

                            D1 = mesh.vol[i,j,k]/AW[phi.gIndex[i,j,k],phi.gIndex[i,j,k]]
                            D2 = mesh.vol[i,j,k-1]/AW[phi.gIndex[i,j,k-1],phi.gIndex[i,j,k-1]]
                            Db = general_interpolation(
                                [mesh.dz[k]; mesh.dz[k-1]], [D1; D2];
                                interpolation = interpolation
                            )

                            num = rho * Db * (mesh.dx[i] * mesh.dy[j])
                            den = mesh.z[k] - mesh.z[k-1]
                            ab = num / den

                            b5 = rho * velocityW[i,j,k] * (mesh.dx[i] * mesh.dy[j])

                            A[id, phi.gIndex[i,j,k-1]] = (-1.0) * ab
                            b[id] += b5
                        else
                            b5 = material.ρ[i,j,k] * velocityW[i,j,k] * (mesh.dx[i] * mesh.dy[j])

                            b[id] += b5
                        end

                        #Top Coefficents
                        if (k != mesh.n1) && (phi.onoff[i,j,k+1])
                            rho = density_interpolation(
                                [mesh.dz[k]; mesh.dz[k+1]], [material.ρ[i,j,k]; material.ρ[i,j,k+1]];
                                interpolation = interpolation
                            )
                            D1 = mesh.vol[i,j,k]/AW[phi.gIndex[i,j,k],phi.gIndex[i,j,k]]
                            D2 = mesh.vol[i,j,k+1]/AW[phi.gIndex[i,j,k+1],phi.gIndex[i,j,k+1]]
                            Dt = general_interpolation(
                                [mesh.dz[k]; mesh.dz[k+1]], [D1; D2];
                                interpolation = interpolation
                            )

                            num = rho * Dt * (mesh.dx[i] * mesh.dy[j])
                            den = mesh.z[k+1] - mesh.z[k]
                            at = num / den

                            b6 = -1.0 * rho * velocityW[i,j,k+1] * (mesh.dx[i] * mesh.dy[j])

                            A[id, phi.gIndex[i,j,k+1]] = (-1.0) * at
                            b[id] += b6
                        else
                            b6 = -1.0 * material.ρ[i,j,k] * velocityW[i,j,k+1] * (mesh.dx[i] * mesh.dy[j])

                            b[id] += b6
                        end

                        #Center Coefficent
                        A[id, id] = aw + ae + as + an + ab + at

                    end
                end
            end
        end

    end

    if sparse
        A = dropzeros(A)
    end

    return A, b
end

function discretize_SIMPLE_PressureCorrection(
    AU::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,2},
    },
    AV::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,2},
    },
    AW::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,2},
    },
    vel::CSVelocity3D,
    phi::CSPhi3D,
    material::UnionCSConstantMaterial,
    mesh::UnionCSMesh3D;
    velocityU::Array{<:AbstractFloat,3} = vel.fValues.uFace,
    velocityV::Array{<:AbstractFloat,3} = vel.fValues.vFace,
    velocityW::Array{<:AbstractFloat,3} = vel.fValues.wFace,
    T::Type{<:AbstractFloat} = Float64,
    threads::Bool = false,
    sparse::Bool = true,
    scheme::Signed = 2,
    interpolation::Signed = 1,
)
    n_equations = maximum_globalIndex(phi)

    if sparse
        A = spzeros(T, n_equations, n_equations)
    elseif !sparse
        A = zeros(T, n_equations, n_equations)
    end

    b = zeros(T, n_equations)

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            for j in 1:mesh.m1
                for k in 1:mesh.n1
                    if phi.onoff[i,j,k]
                        id = phi.gIndex[i,j,k]

                        #Auxiliar variables
                        ac = 0.0
                        b0 = 0.0
                        aw = 0.0
                        Dw = 0.0
                        b1 = 0.0
                        ae = 0.0
                        De = 0.0
                        b2 = 0.0
                        as = 0.0
                        Ds = 0.0
                        b3 = 0.0
                        an = 0.0
                        Dn = 0.0
                        b4 = 0.0
                        ab = 0.0
                        Db = 0.0
                        b5 = 0.0
                        at = 0.0
                        Dt = 0.0
                        b6 = 0.0

                        #West Coefficents
                        if (i != 1) && (phi.onoff[i-1,j,k])
                            rho = material.ρ

                            D1 = mesh.vol[i,j,k]/AU[phi.gIndex[i,j,k],phi.gIndex[i,j,k]]
                            D2 = mesh.vol[i-1,j,k]/AU[phi.gIndex[i-1,j,k],phi.gIndex[i-1,j,k]]
                            Dw = general_interpolation(
                                [mesh.dx[i]; mesh.dx[i-1]], [D1; D2];
                                interpolation = interpolation
                            )

                            num = rho * Dw * (mesh.dy[j] * mesh.dz[k])
                            den = mesh.x[i] - mesh.x[i-1]
                            aw = num / den

                            b1 = rho * velocityU[i,j,k] * (mesh.dy[j] * mesh.dz[k])

                            A[id, phi.gIndex[i-1,j,k]] = (-1.0) * aw
                            b[id] += b1
                        else
                            b1 = material.ρ * velocityU[i,j,k] * (mesh.dy[j] * mesh.dy[j])

                            b[id] += b1
                        end

                        #East Coefficents
                        if (i != mesh.l1)  && (phi.onoff[i+1,j,k])
                            rho = material.ρ

                            D1 = mesh.vol[i,j,k]/AU[phi.gIndex[i,j,k],phi.gIndex[i,j,k]]
                            D2 = mesh.vol[i+1,j,k]/AU[phi.gIndex[i+1,j,k],phi.gIndex[i+1,j,k]]
                            De = general_interpolation(
                                [mesh.dx[i]; mesh.dx[i+1]], [D1; D2];
                                interpolation = interpolation
                            )

                            num = rho * De * (mesh.dy[j] * mesh.dz[k])
                            den = mesh.x[i+1] - mesh.x[i]
                            ae = num / den

                            b2 = -1.0 * rho * velocityU[i+1,j,k] * (mesh.dy[j] * mesh.dz[k])

                            A[id, phi.gIndex[i+1,j,k]] = (-1.0) * ae
                            b[id] += b2
                        else
                            b2 = -1.0 * material.ρ * velocityU[i+1,j,k] * (mesh.dy[j] * mesh.dy[j])

                            b[id] += b2
                        end

                        #South Coefficents
                        if (j != 1) && (phi.onoff[i,j-1,k])
                            rho = material.ρ

                            D1 = mesh.vol[i,j,k]/AV[phi.gIndex[i,j,k],phi.gIndex[i,j,k]]
                            D2 = mesh.vol[i,j-1,k]/AV[phi.gIndex[i,j-1,k],phi.gIndex[i,j-1,k]]
                            Ds = general_interpolation(
                                [mesh.dy[j]; mesh.dy[j-1]], [D1; D2];
                                interpolation = interpolation
                            )

                            num = rho * Ds * (mesh.dx[i] * mesh.dz[k])
                            den = mesh.y[j] - mesh.y[j-1]
                            as = num / den

                            b3 = rho * velocityV[i,j,k] * (mesh.dx[i] * mesh.dz[k])

                            A[id, phi.gIndex[i,j-1,k]] = (-1.0) * as
                            b[id] += b3
                        else
                            b3 = material.ρ * velocityV[i,j,k] * (mesh.dx[i] * mesh.dz[k])

                            b[id] += b3
                        end

                        #North Coefficents
                        if (j != mesh.l1) && (phi.onoff[i,j+1,k])
                            rho = material.ρ

                            D1 = mesh.vol[i,j,k]/AV[phi.gIndex[i,j,k],phi.gIndex[i,j,k]]
                            D2 = mesh.vol[i,j+1,k]/AV[phi.gIndex[i,j+1,k],phi.gIndex[i,j+1,k]]
                            Dn = general_interpolation(
                                [mesh.dy[j]; mesh.dy[j+1]], [D1; D2];
                                interpolation = interpolation
                            )

                            num = rho * Dn * (mesh.dx[i] * mesh.dz[k])
                            den = mesh.y[j+1] - mesh.y[j]
                            an = num / den

                            b4 = -1.0 * rho * velocityV[i,j+1,k] * (mesh.dx[i] * mesh.dz[k])

                            A[id, phi.gIndex[i,j+1,k]] = (-1.0) * an
                            b[id] += b4
                        else
                            b4 = -1.0 * material.ρ * velocityV[i,j+1,k] * (mesh.dx[i] * mesh.dz[k])

                            b[id] += b4
                        end

                        #Bottom Coefficents
                        if (k != 1) && (phi.onoff[i,j,k-1])
                            rho = material.ρ

                            D1 = mesh.vol[i,j,k]/AW[phi.gIndex[i,j,k],phi.gIndex[i,j,k]]
                            D2 = mesh.vol[i,j,k-1]/AW[phi.gIndex[i,j,k-1],phi.gIndex[i,j,k-1]]
                            Db = general_interpolation(
                                [mesh.dz[k]; mesh.dz[k-1]], [D1; D2];
                                interpolation = interpolation
                            )

                            num = rho * Db * (mesh.dx[i] * mesh.dy[j])
                            den = mesh.z[k] - mesh.z[k-1]
                            ab = num / den

                            b5 = rho * velocityW[i,j,k] * (mesh.dx[i] * mesh.dy[j])

                            A[id, phi.gIndex[i,j,k-1]] = (-1.0) * ab
                            b[id] += b5
                        else
                            b5 = material.ρ * velocityW[i,j,k] * (mesh.dx[i] * mesh.dy[j])

                            b[id] += b5
                        end

                        #Top Coefficents
                        if (k != mesh.n1) && (phi.onoff[i,j,k+1])
                            rho = material.ρ

                            D1 = mesh.vol[i,j,k]/AW[phi.gIndex[i,j,k],phi.gIndex[i,j,k]]
                            D2 = mesh.vol[i,j,k+1]/AW[phi.gIndex[i,j,k+1],phi.gIndex[i,j,k+1]]
                            Dt = general_interpolation(
                                [mesh.dz[k]; mesh.dz[k+1]], [D1; D2];
                                interpolation = interpolation
                            )

                            num = rho * Dt * (mesh.dx[i] * mesh.dy[j])
                            den = mesh.z[k+1] - mesh.z[k]
                            at = num / den

                            b6 = -1.0 * rho * velocityW[i,j,k+1] * (mesh.dx[i] * mesh.dy[j])

                            A[id, phi.gIndex[i,j,k+1]] = (-1.0) * at
                            b[id] += b6
                        else
                            b6 = -1.0 * material.ρ * velocityW[i,j,k+1] * (mesh.dx[i] * mesh.dy[j])

                            b[id] += b6
                        end

                        #Center Coefficent
                        A[id, id] = aw + ae + as + an + ab + at

                    end
                end
            end
        end

    elseif !threads
        for i in 1:mesh.l1
            for j in 1:mesh.m1
                for k in 1:mesh.n1
                    if phi.onoff[i,j,k]
                        id = phi.gIndex[i,j,k]

                        #Auxiliar variables
                        ac = 0.0
                        b0 = 0.0
                        aw = 0.0
                        Dw = 0.0
                        b1 = 0.0
                        ae = 0.0
                        De = 0.0
                        b2 = 0.0
                        as = 0.0
                        Ds = 0.0
                        b3 = 0.0
                        an = 0.0
                        Dn = 0.0
                        b4 = 0.0
                        ab = 0.0
                        Db = 0.0
                        b5 = 0.0
                        at = 0.0
                        Dt = 0.0
                        b6 = 0.0

                        #West Coefficents
                        if (i != 1) && (phi.onoff[i-1,j,k])
                            rho = material.ρ

                            D1 = mesh.vol[i,j,k]/AU[phi.gIndex[i,j,k],phi.gIndex[i,j,k]]
                            D2 = mesh.vol[i-1,j,k]/AU[phi.gIndex[i-1,j,k],phi.gIndex[i-1,j,k]]
                            Dw = general_interpolation(
                                [mesh.dx[i]; mesh.dx[i-1]], [D1; D2];
                                interpolation = interpolation
                            )

                            num = rho * Dw * (mesh.dy[j] * mesh.dz[k])
                            den = mesh.x[i] - mesh.x[i-1]
                            aw = num / den

                            b1 = rho * velocityU[i,j,k] * (mesh.dy[j] * mesh.dz[k])

                            A[id, phi.gIndex[i-1,j,k]] = (-1.0) * aw
                            b[id] += b1
                        else
                            b1 = material.ρ * velocityU[i,j,k] * (mesh.dy[j] * mesh.dy[j])

                            b[id] += b1
                        end

                        #East Coefficents
                        if (i != mesh.l1)  && (phi.onoff[i+1,j,k])
                            rho = material.ρ

                            D1 = mesh.vol[i,j,k]/AU[phi.gIndex[i,j,k],phi.gIndex[i,j,k]]
                            D2 = mesh.vol[i+1,j,k]/AU[phi.gIndex[i+1,j,k],phi.gIndex[i+1,j,k]]
                            De = general_interpolation(
                                [mesh.dx[i]; mesh.dx[i+1]], [D1; D2];
                                interpolation = interpolation
                            )

                            num = rho * De * (mesh.dy[j] * mesh.dz[k])
                            den = mesh.x[i+1] - mesh.x[i]
                            ae = num / den

                            b2 = -1.0 * rho * velocityU[i+1,j,k] * (mesh.dy[j] * mesh.dz[k])

                            A[id, phi.gIndex[i+1,j,k]] = (-1.0) * ae
                            b[id] += b2
                        else
                            b2 = -1.0 * material.ρ * velocityU[i+1,j,k] * (mesh.dy[j] * mesh.dy[j])

                            b[id] += b2
                        end

                        #South Coefficents
                        if (j != 1) && (phi.onoff[i,j-1,k])
                            rho = material.ρ

                            D1 = mesh.vol[i,j,k]/AV[phi.gIndex[i,j,k],phi.gIndex[i,j,k]]
                            D2 = mesh.vol[i,j-1,k]/AV[phi.gIndex[i,j-1,k],phi.gIndex[i,j-1,k]]
                            Ds = general_interpolation(
                                [mesh.dy[j]; mesh.dy[j-1]], [D1; D2];
                                interpolation = interpolation
                            )

                            num = rho * Ds * (mesh.dx[i] * mesh.dz[k])
                            den = mesh.y[j] - mesh.y[j-1]
                            as = num / den

                            b3 = rho * velocityV[i,j,k] * (mesh.dx[i] * mesh.dz[k])

                            A[id, phi.gIndex[i,j-1,k]] = (-1.0) * as
                            b[id] += b3
                        else
                            b3 = material.ρ * velocityV[i,j,k] * (mesh.dx[i] * mesh.dz[k])

                            b[id] += b3
                        end

                        #North Coefficents
                        if (j != mesh.l1) && (phi.onoff[i,j+1,k])
                            rho = material.ρ

                            D1 = mesh.vol[i,j,k]/AV[phi.gIndex[i,j,k],phi.gIndex[i,j,k]]
                            D2 = mesh.vol[i,j+1,k]/AV[phi.gIndex[i,j+1,k],phi.gIndex[i,j+1,k]]
                            Dn = general_interpolation(
                                [mesh.dy[j]; mesh.dy[j+1]], [D1; D2];
                                interpolation = interpolation
                            )

                            num = rho * Dn * (mesh.dx[i] * mesh.dz[k])
                            den = mesh.y[j+1] - mesh.y[j]
                            an = num / den

                            b4 = -1.0 * rho * velocityV[i,j+1,k] * (mesh.dx[i] * mesh.dz[k])

                            A[id, phi.gIndex[i,j+1,k]] = (-1.0) * an
                            b[id] += b4
                        else
                            b4 = -1.0 * material.ρ * velocityV[i,j+1,k] * (mesh.dx[i] * mesh.dz[k])

                            b[id] += b4
                        end

                        #Bottom Coefficents
                        if (k != 1) && (phi.onoff[i,j,k-1])
                            rho = material.ρ

                            D1 = mesh.vol[i,j,k]/AW[phi.gIndex[i,j,k],phi.gIndex[i,j,k]]
                            D2 = mesh.vol[i,j,k-1]/AW[phi.gIndex[i,j,k-1],phi.gIndex[i,j,k-1]]
                            Db = general_interpolation(
                                [mesh.dz[k]; mesh.dz[k-1]], [D1; D2];
                                interpolation = interpolation
                            )

                            num = rho * Db * (mesh.dx[i] * mesh.dy[j])
                            den = mesh.z[k] - mesh.z[k-1]
                            ab = num / den

                            b5 = rho * velocityW[i,j,k] * (mesh.dx[i] * mesh.dy[j])

                            A[id, phi.gIndex[i,j,k-1]] = (-1.0) * ab
                            b[id] += b5
                        else
                            b5 = material.ρ * velocityW[i,j,k] * (mesh.dx[i] * mesh.dy[j])

                            b[id] += b5
                        end

                        #Top Coefficents
                        if (k != mesh.n1) && (phi.onoff[i,j,k+1])
                            rho = material.ρ

                            D1 = mesh.vol[i,j,k]/AW[phi.gIndex[i,j,k],phi.gIndex[i,j,k]]
                            D2 = mesh.vol[i,j,k+1]/AW[phi.gIndex[i,j,k+1],phi.gIndex[i,j,k+1]]
                            Dt = general_interpolation(
                                [mesh.dz[k]; mesh.dz[k+1]], [D1; D2];
                                interpolation = interpolation
                            )

                            num = rho * Dt * (mesh.dx[i] * mesh.dy[j])
                            den = mesh.z[k+1] - mesh.z[k]
                            at = num / den

                            b6 = -1.0 * rho * velocityW[i,j,k+1] * (mesh.dx[i] * mesh.dy[j])

                            A[id, phi.gIndex[i,j,k+1]] = (-1.0) * at
                            b[id] += b6
                        else
                            b6 = -1.0 * material.ρ * velocityW[i,j,k+1] * (mesh.dx[i] * mesh.dy[j])

                            b[id] += b6
                        end

                        #Center Coefficent
                        A[id, id] = aw + ae + as + an + ab + at

                    end
                end
            end
        end

    end

    if sparse
        A = dropzeros(A)
    end

    return A, b
end
