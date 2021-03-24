"""

"""
function _SIMPLE_PressureCorrection_Coefficients_ end

function _SIMPLE_PressureCorrection_Coefficients_(
    AU::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,2},
    },
    i::Signed,
    velocityU::Array{<:AbstractFloat,1},
    phi::CSPhi1D,
    material::CSMaterial1D,
    mesh::UnionCSMesh1D;
    T::Type{<:AbstractFloat} = Float64,
    sparse::Bool = true,
    interpolation::Signed = 1,
)
    n_equations = maximum_globalIndex(phi)

    if sparse
        A = spzeros(T, 1, n_equations)
    elseif !sparse
        A = zeros(T, 1, n_equations)
    end

    b = 0.0

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

        A[1, phi.gIndex[i-1]] = (-1.0) * aw
        b += b1

    else
        b1 = material.ρ[i] * velocityU[i] * 1.0

        b += b1
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

        A[1, phi.gIndex[i+1]] = (-1.0) * ae
        b += b2

    else
        b2 = -1.0 * material.ρ[i] * velocityU[i+1] * 1.0

        b += b2
    end

    #Center Coefficent
    A[1, phi.gIndex[i]] = ae + aw

    return A, b
end

function _SIMPLE_PressureCorrection_Coefficients_(
    AU::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,2},
    },
    i::Signed,
    velocityU::Array{<:AbstractFloat,1},
    phi::CSPhi1D,
    material::UnionCSConstantMaterial,
    mesh::UnionCSMesh1D;
    T::Type{<:AbstractFloat} = Float64,
    sparse::Bool = true,
    interpolation::Signed = 1,
)
    n_equations = maximum_globalIndex(phi)

    if sparse
        A = spzeros(T, 1, n_equations)
    elseif !sparse
        A = zeros(T, 1, n_equations)
    end

    b = 0.0

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

        A[1, phi.gIndex[i-1]] = (-1.0) * aw
        b += b1
    else
        b1 = material.ρ * velocityU[i] * 1.0

        b += b1
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

        A[1, phi.gIndex[i+1]] = (-1.0) * ae
        b += b2
    else
        b2 = -1.0 * material.ρ * velocityU[i+1] * 1.0

        b += b2
    end

    #Center Coefficent
    A[1, phi.gIndex[i]] = ae + aw

    return A, b
end

function _SIMPLE_PressureCorrection_Coefficients_(
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
    i::Signed,
    j::Signed,
    velocityU::Array{<:AbstractFloat,2},
    velocityV::Array{<:AbstractFloat,2},
    phi::CSPhi2D,
    material::CSMaterial2D,
    mesh::UnionCSMesh2D;
    T::Type{<:AbstractFloat} = Float64,
    sparse::Bool = true,
    interpolation::Signed = 2,
)

    n_equations = maximum_globalIndex(phi)

    if sparse
        A = spzeros(T, 1, n_equations)
    elseif !sparse
        A = zeros(T, 1, n_equations)
    end

    b = 0.0

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

        A[1, phi.gIndex[i-1,j]] = (-1.0) * aw
        b += b1
    else
        b1 = material.ρ[i,j] * velocityU[i,j] * mesh.dy[j]

        b += b1
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

        A[1, phi.gIndex[i+1,j]] = (-1.0) * ae
        b += b2
    else
        b2 = -1.0 * material.ρ[i,j] * velocityU[i+1,j] * mesh.dy[j]

        b += b2
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

        A[1, phi.gIndex[i,j-1]] = (-1.0) * as
        b += b3
    else
        b3 = material.ρ[i,j] * velocityV[i,j] * mesh.dx[i]

        b += b3
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

        A[1, phi.gIndex[i,j+1]] = (-1.0) * an
        b += b4
    else
        b4 = -1.0 * material.ρ[i,j] * velocityV[i,j+1] * mesh.dx[i]

        b += b4
    end

    #Center Coefficent
    A[1, phi.gIndex[i,j]] = ae + aw + as + an

    return A, b
end

function _SIMPLE_PressureCorrection_Coefficients_(
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
    i::Signed,
    j::Signed,
    velocityU::Array{<:AbstractFloat,2},
    velocityV::Array{<:AbstractFloat,2},
    phi::CSPhi2D,
    material::UnionCSConstantMaterial,
    mesh::UnionCSMesh2D;
    T::Type{<:AbstractFloat} = Float64,
    sparse::Bool = true,
    interpolation::Signed = 2,
)

    n_equations = maximum_globalIndex(phi)

    if sparse
        A = spzeros(T, 1, n_equations)
    elseif !sparse
        A = zeros(T, 1, n_equations)
    end

    b = 0.0

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

        A[1, phi.gIndex[i-1,j]] = (-1.0) * aw
        b += b1
    else
        b1 = material.ρ * velocityU[i,j] * mesh.dy[j]

        b += b1
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

        A[1, phi.gIndex[i+1,j]] = (-1.0) * ae
        b += b2
    else
        b2 = -1.0 * material.ρ * velocityU[i+1,j] * mesh.dy[j]

        b += b2
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

        A[1, phi.gIndex[i,j-1]] = (-1.0) * as
        b += b3
    else
        b3 = material.ρ * velocityV[i,j] * mesh.dx[i]

        b += b3
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

        A[1, phi.gIndex[i,j+1]] = (-1.0) * an
        b += b4
    else
        b4 = -1.0 * material.ρ * velocityV[i,j+1] * mesh.dx[i]

        b += b4
    end

    #Center Coefficent
    A[1, phi.gIndex[i,j]] = ae + aw + as + an

    return A, b
end

function _SIMPLE_PressureCorrection_Coefficients_(
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
    i::Signed,
    j::Signed,
    k::Signed,
    velocityU::Array{<:AbstractFloat,3},
    velocityV::Array{<:AbstractFloat,3},
    velocityW::Array{<:AbstractFloat,3},
    phi::CSPhi3D,
    material::CSMaterial3D,
    mesh::UnionCSMesh3D;
    T::Type{<:AbstractFloat} = Float64,
    sparse::Bool = true,
    interpolation::Signed = 2,
)
    n_equations = maximum_globalIndex(phi)

    if sparse
        A = spzeros(T, 1, n_equations)
    elseif !sparse
        A = zeros(T, 1, n_equations)
    end

    b = 0.0

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

        A[1, phi.gIndex[i-1,j,k]] = (-1.0) * aw
        b += b1
    else
        b1 = material.ρ[i,j,k] * velocityU[i,j,k] * (mesh.dy[j] * mesh.dz[k])

        b += b1
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

        A[1, phi.gIndex[i+1,j,k]] = (-1.0) * ae
        b += b2
    else
        b2 = -1.0 * material.ρ[i,j,k] * velocityU[i+1,j,k] * (mesh.dy[j] * mesh.dz[k])

        b += b2
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

        A[1, phi.gIndex[i,j-1,k]] = (-1.0) * as
        b += b3
    else
        b3 = material.ρ[i,j,k] * velocityV[i,j,k] * (mesh.dx[i] * mesh.dz[k])

        b += b3
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

        A[1, phi.gIndex[i,j+1,k]] = (-1.0) * an
        b += b4
    else
        b4 = -1.0 * material.ρ[i,j,k] * velocityV[i,j+1,k] * (mesh.dx[i] * mesh.dz[k])

        b += b4
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

        A[1, phi.gIndex[i,j,k-1]] = (-1.0) * ab
        b += b5
    else
        b5 = material.ρ[i,j,k] * velocityW[i,j,k] * (mesh.dx[i] * mesh.dy[j])

        b += b5
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

        A[1, phi.gIndex[i,j,k+1]] = (-1.0) * at
        b += b6
    else
        b6 = -1.0 * material.ρ[i,j,k] * velocityW[i,j,k+1] * (mesh.dx[i] * mesh.dy[j])

        b += b6
    end

    #Center Coefficent
    A[1, phi.gIndex[i,j,k]] = aw + ae + as + an + ab + at

    return A, b
end

function _SIMPLE_PressureCorrection_Coefficients_(
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
    i::Signed,
    j::Signed,
    k::Signed,
    velocityU::Array{<:AbstractFloat,3},
    velocityV::Array{<:AbstractFloat,3},
    velocityW::Array{<:AbstractFloat,3},
    phi::CSPhi3D,
    material::UnionCSConstantMaterial,
    mesh::UnionCSMesh3D;
    T::Type{<:AbstractFloat} = Float64,
    sparse::Bool = true,
    interpolation::Signed = 2,
)
    n_equations = maximum_globalIndex(phi)

    if sparse
        A = spzeros(T, 1, n_equations)
    elseif !sparse
        A = zeros(T, 1, n_equations)
    end

    b = 0.0

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

        A[1, phi.gIndex[i-1,j,k]] = (-1.0) * aw
        b += b1
    else
        b1 = material.ρ * velocityU[i,j,k] * (mesh.dy[j] * mesh.dy[j])

        b += b1
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

        A[1, phi.gIndex[i+1,j,k]] = (-1.0) * ae
        b += b2
    else
        b2 = -1.0 * material.ρ * velocityU[i+1,j,k] * (mesh.dy[j] * mesh.dy[j])

        b += b2
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

        A[1, phi.gIndex[i,j-1,k]] = (-1.0) * as
        b += b3
    else
        b3 = material.ρ * velocityV[i,j,k] * (mesh.dx[i] * mesh.dz[k])

        b += b3
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

        A[1, phi.gIndex[i,j+1,k]] = (-1.0) * an
        b += b4
    else
        b4 = -1.0 * material.ρ * velocityV[i,j+1,k] * (mesh.dx[i] * mesh.dz[k])

        b += b4
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

        A[1, phi.gIndex[i,j,k-1]] = (-1.0) * ab
        b += b5
    else
        b5 = material.ρ * velocityW[i,j,k] * (mesh.dx[i] * mesh.dy[j])

        b += b5
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

        A[1, phi.gIndex[i,j,k+1]] = (-1.0) * at
        b += b6
    else
        b6 = -1.0 * material.ρ * velocityW[i,j,k+1] * (mesh.dx[i] * mesh.dy[j])

        b += b6
    end

    #Center Coefficent
    A[1, phi.gIndex[i,j,k]] = aw + ae + as + an + ab + at

    return A, b
end
