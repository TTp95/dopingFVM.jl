"""

"""
function _discretize_convection_centralDifference_ end

function _discretize_convection_centralDifference_(
    i::Signed,
    velocityU::Array{AbstractFloat,1},
    phi::CSPhi1D,
    bounds::Dict{String,BoundsStructured},
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
    awc = 0.0
    b1 = 0.0
    ae = 0.0
    aec = 0.0
    b2 = 0.0

    #West Coefficents
    if (i != 1) && (phi.onoff[i-1])
        rho = gamma_interpolation(
            [mesh.dx[i]; mesh.dx[i-1]], [material.ρ[i]; material.ρ[i-1]];
            interpolation = interpolation
        )
        @inbounds aw, awc, b1 = _convection_centralDifference_neighbors_(
            rho, velocityU[i], [mesh.dx[i]; mesh.dx[i-1]], (1.0)
        )
        A[1, phi.gIndex[i-1]] = (-1.0) * aw
        b += b1
    end

    #East Coefficents
    if (i != mesh.l1)  && (phi.onoff[i+1])
        rho = gamma_interpolation(
            [mesh.dx[i]; mesh.dx[i+1]], [material.ρ[i]; material.ρ[i+1]];
            interpolation = interpolation
        )
        @inbounds ae, aec, b2 = _convection_centralDifference_neighbors_(
        rho, velocityU[i+1], [mesh.dx[i]; mesh.dx[i+1]], (1.0)
        )
        A[1, phi.gIndex[i+1]] = ae
        b += b2
    end

    #Center Coefficent
    if phi.bounds[i]
        @inbounds ac, b0 = _diffusion_centralDifference_central_(
            i,
            velocityU,
            phi,
            bounds,
            mesh;
            T = T,
        )
        A[1, phi.gIndex[i]] = ac + aec - awc
        b += b0
    else
        A[1, phi.gIndex[i]] = ac + aec - awc
    end

    return A, b
end

function _discretize_diffusion_centralDifference_(
    i::Signed,
    velocityU::Array{AbstractFloat,1},
    phi::CSPhi1D,
    bounds::Dict{String,BoundsStructured},
    material::UnionCSConstantMaterial,
    mesh::UnionCSMesh1D;
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
    awc = 0.0
    b1 = 0.0
    ae = 0.0
    aec = 0.0
    b2 = 0.0

    #West Coefficents
    if (i != 1) && (phi.onoff[i-1])
        rho = material.ρ
        @inbounds aw, awc, b1 = _convection_centralDifference_neighbors_(
            rho, velocityU[i], [mesh.dx[i]; mesh.dx[i-1]], (1.0)
        )
        A[1, phi.gIndex[i-1]] = (-1.0) * aw
        b += b1
    end

    #East Coefficents
    if (i != mesh.l1)  && (phi.onoff[i+1])
        rho = material.ρ
        @inbounds ae, aec, b2 = _convection_centralDifference_neighbors_(
        rho, velocityU[i+1], [mesh.dx[i]; mesh.dx[i+1]], (1.0)
        )
        A[1, phi.gIndex[i+1]] = ae
        b += b2
    end

    #Center Coefficent
    if phi.bounds[i]
        @inbounds ac, b0 = _diffusion_centralDifference_central_(
            i,
            velocityU,
            phi,
            bounds,
            mesh;
            T = T,
        )
        A[1, phi.gIndex[i]] = ac + aec - awc
        b += b0
    else
        A[1, phi.gIndex[i]] = ac + aec - awc
    end

    return A, b
end

function _discretize_diffusion_centralDifference_(
    i::Signed,
    j::Signed,
    velocityU::Array{AbstractFloat,2},
    velocityV::Array{AbstractFloat,2},
    phi::CSPhi2D,
    bounds::Dict{String,BoundsStructured},
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
    awc = 0.0
    b1 = 0.0
    ae = 0.0
    aec = 0.0
    b2 = 0.0
    as = 0.0
    asc = 0.0
    b3 = 0.0
    an = 0.0
    anc = 0.0
    b4 = 0.0

    #West Coefficents
    if (i != 1) && (phi.onoff[i-1,j])
        rho = gamma_interpolation(
            [mesh.dx[i]; mesh.dx[i-1]], [material.ρ[i,j]; material.ρ[i-1,j]];
            interpolation = interpolation
        )
        @inbounds aw, awc, b1 = _convection_centralDifference_neighbors_(
            rho, velocityU[i,j], [mesh.dx[i]; mesh.dx[i-1]], (mesh.dy[j])
        )
        A[1, phi.gIndex[i-1,j]] = (-1.0) * aw
        b += b1
    end

    #East Coefficents
    if (i != mesh.l1)  && (phi.onoff[i+1,j])
        rho = gamma_interpolation(
            [mesh.dx[i]; mesh.dx[i+1]], [material.ρ[i,j]; material.ρ[i+1,j]];
            interpolation = interpolation
        )
        @inbounds ae, aec, b2 = _convection_centralDifference_neighbors_(
        rho, velocityU[i+1,j], [mesh.dx[i]; mesh.dx[i+1]], (mesh.dy[j])
        )
        A[1, phi.gIndex[i+1,j]] = ae
        b += b2
    end

    #South Coefficents
    if (j != 1) && (phi.onoff[i,j-1])
        rho = gamma_interpolation(
            [mesh.dy[j]; mesh.dy[j-1]], [material.ρ[i,j]; material.ρ[i,j-1]];
            interpolation = interpolation
        )
        @inbounds as, asc, b3 = _convection_centralDifference_neighbors_(
            rho, velocityV[i,j], [mesh.dy[j]; mesh.dy[j-1]], (mesh.dx[i])
        )
        A[1, phi.gIndex[i,j-1]] = (-1.0) * as
        b += b3
    end

    #North Coefficents
    if (j != mesh.m1) && (phi.onoff[i,j+1])
        rho = gamma_interpolation(
            [mesh.dy[j]; mesh.dy[j+1]], [material.ρ[i,j]; material.ρ[i,j+1]];
            interpolation = interpolation
        )
        @inbounds an, anc, b4 = _convection_centralDifference_neighbors_(
            rho, velocityV[i,j+1], [mesh.dy[j]; mesh.dy[j+1]], (mesh.dx[i])
        )
        A[1, phi.gIndex[i,j+1]] = an
        b += b4
    end

    #Center Coefficent
    if phi.bounds[i,j]
        @inbounds ac, b0 = _diffusion_centralDifference_central_(
            i,
            j,
            velocityU,
            velocityV,
            phi,
            bounds,
            mesh;
            T = T,
        )
        A[1, phi.gIndex[i,j]] = ac + aec - awc + anc - asc
        b += b0
    else
        A[1, phi.gIndex[i,j]] = ac + aec - awc + anc - asc
    end

    return A, b
end

function _discretize_diffusion_centralDifference_(
    i::Signed,
    j::Signed,
    velocityU::Array{AbstractFloat,2},
    velocityV::Array{AbstractFloat,2},
    phi::CSPhi2D,
    bounds::Dict{String,BoundsStructured},
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
    awc = 0.0
    b1 = 0.0
    ae = 0.0
    aec = 0.0
    b2 = 0.0
    as = 0.0
    asc = 0.0
    b3 = 0.0
    an = 0.0
    anc = 0.0
    b4 = 0.0

    #West Coefficents
    if (i != 1) && (phi.onoff[i-1,j])
        rho = material.ρ
        @inbounds aw, awc, b1 = _convection_centralDifference_neighbors_(
            rho, velocityU[i,j], [mesh.dx[i]; mesh.dx[i-1]], (mesh.dy[j])
        )
        A[1, phi.gIndex[i-1,j]] = (-1.0) * aw
        b += b1
    end

    #East Coefficents
    if (i != mesh.l1)  && (phi.onoff[i+1,j])
        rho = material.ρ
        @inbounds ae, aec, b2 = _convection_centralDifference_neighbors_(
        rho, velocityU[i+1,j], [mesh.dx[i]; mesh.dx[i+1]], (mesh.dy[j])
        )
        A[1, phi.gIndex[i+1,j]] = ae
        b += b2
    end

    #South Coefficents
    if (j != 1) && (phi.onoff[i,j-1])
        rho = material.ρ
        @inbounds as, asc, b3 = _convection_centralDifference_neighbors_(
            rho, velocityV[i,j], [mesh.dy[j]; mesh.dy[j-1]], (mesh.dx[i])
        )
        A[1, phi.gIndex[i,j-1]] = (-1.0) * as
        b += b3
    end

    #North Coefficents
    if (j != mesh.m1) && (phi.onoff[i,j+1])
        rho = material.ρ
        @inbounds an, anc, b4 = _convection_centralDifference_neighbors_(
            rho, velocityV[i,j+1], [mesh.dy[j]; mesh.dy[j+1]], (mesh.dx[i])
        )
        A[1, phi.gIndex[i,j+1]] = an
        b += b4
    end

    #Center Coefficent
    if phi.bounds[i,j]
        @inbounds ac, b0 = _diffusion_centralDifference_central_(
            i,
            j,
            velocityU,
            velocityV,
            phi,
            bounds,
            mesh;
            T = T,
        )
        A[1, phi.gIndex[i,j]] = ac + aec - awc + anc - asc
        b += b0
    else
        A[1, phi.gIndex[i,j]] = ac + aec - awc + anc - asc
    end

    return A, b
end

function _discretize_diffusion_centralDifference_(
    i::Signed,
    j::Signed,
    k::Signed,
    velocityU::Array{AbstractFloat,3},
    velocityV::Array{AbstractFloat,3},
    velocityW::Array{AbstractFloat,3},
    phi::CSPhi3D,
    bounds::Dict{String,BoundsStructured},
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
    awc = 0.0
    b1 = 0.0
    ae = 0.0
    aec = 0.0
    b2 = 0.0
    as = 0.0
    asc = 0.0
    b3 = 0.0
    an = 0.0
    anc = 0.0
    b4 = 0.0
    ab = 0.0
    abc = 0.0
    b5 = 0.0
    at = 0.0
    atc = 0.0
    b6 = 0.0

    #West Coefficents
    if (i != 1) && (phi.onoff[i-1,j,k])
        rho = gamma_interpolation(
            [mesh.dx[i]; mesh.dx[i-1]], [material.ρ[i,j,k]; material.ρ[i-1,j,k]];
            interpolation = interpolation
        )
        @inbounds aw, awc, b1 = _convection_centralDifference_neighbors_(
            rho, velocityU[i,j,k], [mesh.dx[i]; mesh.dx[i-1]], (mesh.dy[j] * mesh.dz[k])
        )
        A[1, phi.gIndex[i-1,j,k]] = (-1.0) * aw
        b += b1
    end

    #East Coefficents
    if (i != mesh.l1)  && (phi.onoff[i+1,j,k])
        rho = gamma_interpolation(
            [mesh.dx[i]; mesh.dx[i+1]], [material.ρ[i,j,k]; material.ρ[i+1,j,k]];
            interpolation = interpolation
        )
        @inbounds ae, aec, b2 = _convection_centralDifference_neighbors_(
        rho, velocityU[i+1,j,k], [mesh.dx[i]; mesh.dx[i+1]], (mesh.dy[j] * mesh.dz[k])
        )
        A[1, phi.gIndex[i+1,j,k]] = ae
        b += b2
    end

    #South Coefficents
    if (j != 1) && (phi.onoff[i,j-1,k])
        rho = gamma_interpolation(
            [mesh.dy[j]; mesh.dy[j-1]], [material.ρ[i,j,k]; material.ρ[i,j-1,k]];
            interpolation = interpolation
        )
        @inbounds as, asc, b3 = _convection_centralDifference_neighbors_(
            rho, velocityV[i,j,k], [mesh.dy[j]; mesh.dy[j-1]], (mesh.dx[i] * mesh.dz[k])
        )
        A[1, phi.gIndex[i,j-1,k]] = (-1.0) * as
        b += b3
    end

    #North Coefficents
    if (j != mesh.l1) && (phi.onoff[i,j+1,k])
        rho = gamma_interpolation(
            [mesh.dy[j]; mesh.dy[j+1]], [material.ρ[i,j,k]; material.ρ[i,j+1,k]];
            interpolation = interpolation
        )
        @inbounds an, anc, b4 = _convection_centralDifference_neighbors_(
            rho, velocityV[i,j+1,k], [mesh.dy[j]; mesh.dy[j+1]], (mesh.dx[i] * mesh.dz[k])
        )
        A[1, phi.gIndex[i,j+1,k]] = an
        b += b4
    end

    #Bottom Coefficents
    if (k != 1) && (phi.onoff[i,j,k-1])
        rho = gamma_interpolation(
            [mesh.dz[k]; mesh.dzk[k-1]], [material.ρ[i,j,k]; material.ρ[i,j,k-1]];
            interpolation = interpolation
        )
        @inbounds ab, abc, b5 = _convection_centralDifference_neighbors_(
            rho, velocityW[i,j,k], [mesh.dz[k]; mesh.dz[k-1]], (mesh.dx[i] * mesh.dy[j])
        )
        A[1, phi.gIndex[i,j,k-1]] = (-1.0) * ab
        b += b5
    end

    #Top Coefficents
    if (k != mesh.n1) && (phi.onoff[i,j,k+1])
        rho = gamma_interpolation(
            [mesh.dz[k]; mesh.dzk[k+1]], [material.ρ[i,j,k]; material.ρ[i,j,k+1]];
            interpolation = interpolation
        )
        @inbounds at, atc, b6 = _convection_centralDifference_neighbors_(
            rho, velocityW[i,j,k+1], [mesh.dz[k]; mesh.dz[k+1]], (mesh.dx[i] * mesh.dy[j])
        )
        A[1, phi.gIndex[i,j,k+1]] = at
        b += b6
    end

    #Center Coefficent
    if phi.bounds[i,j,k]
        @inbounds ac, b0 = _diffusion_centralDifference_central_(
            i,
            j,
            k,
            velocityU,
            velocityV,
            velocityW,
            phi,
            bounds,
            mesh;
            T = T,
        )
        A[1, phi.gIndex[i,j,k]] = ac + aec - awc + anc - asc + atc - abc
        b += b0
    else
        A[1, phi.gIndex[i,j,k]] = ac + aec - awc + anc - asc + atc - abc
    end

    return A, b
end

function _discretize_diffusion_centralDifference_(
    i::Signed,
    j::Signed,
    k::Signed,
    velocityU::Array{AbstractFloat,3},
    velocityV::Array{AbstractFloat,3},
    velocityW::Array{AbstractFloat,3},
    phi::CSPhi3D,
    bounds::Dict{String,BoundsStructured},
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
    awc = 0.0
    b1 = 0.0
    ae = 0.0
    aec = 0.0
    b2 = 0.0
    as = 0.0
    asc = 0.0
    b3 = 0.0
    an = 0.0
    anc = 0.0
    b4 = 0.0
    ab = 0.0
    abc = 0.0
    b5 = 0.0
    at = 0.0
    atc = 0.0
    b6 = 0.0

    #West Coefficents
    if (i != 1) && (phi.onoff[i-1,j,k])
        rho = material.ρ
        @inbounds aw, awc, b1 = _convection_centralDifference_neighbors_(
            rho, velocityU[i,j,k], [mesh.dx[i]; mesh.dx[i-1]], (mesh.dy[j] * mesh.dz[k])
        )
        A[1, phi.gIndex[i-1,j,k]] = (-1.0) * aw
        b += b1
    end

    #East Coefficents
    if (i != mesh.l1)  && (phi.onoff[i+1,j,k])
        rho = material.ρ
        @inbounds ae, aec, b2 = _convection_centralDifference_neighbors_(
        rho, velocityU[i+1,j,k], [mesh.dx[i]; mesh.dx[i+1]], (mesh.dy[j] * mesh.dz[k])
        )
        A[1, phi.gIndex[i+1,j,k]] = ae
        b += b2
    end

    #South Coefficents
    if (j != 1) && (phi.onoff[i,j-1,k])
        rho = material.ρ
        @inbounds as, asc, b3 = _convection_centralDifference_neighbors_(
            rho, velocityV[i,j,k], [mesh.dy[j]; mesh.dy[j-1]], (mesh.dx[i] * mesh.dz[k])
        )
        A[1, phi.gIndex[i,j-1,k]] = (-1.0) * as
        b += b3
    end

    #North Coefficents
    if (j != mesh.l1) && (phi.onoff[i,j+1,k])
        rho = material.ρ
        @inbounds an, anc, b4 = _convection_centralDifference_neighbors_(
            rho, velocityV[i,j+1,k], [mesh.dy[j]; mesh.dy[j+1]], (mesh.dx[i] * mesh.dz[k])
        )
        A[1, phi.gIndex[i,j+1,k]] = an
        b += b4
    end

    #Bottom Coefficents
    if (k != 1) && (phi.onoff[i,j,k-1])
        rho = material.ρ
        @inbounds ab, abc, b5 = _convection_centralDifference_neighbors_(
            rho, velocityW[i,j,k], [mesh.dz[k]; mesh.dz[k-1]], (mesh.dx[i] * mesh.dy[j])
        )
        A[1, phi.gIndex[i,j,k-1]] = (-1.0) * ab
        b += b5
    end

    #Top Coefficents
    if (k != mesh.n1) && (phi.onoff[i,j,k+1])
        rho = material.ρ
        @inbounds at, atc, b6 = _convection_centralDifference_neighbors_(
            rho, velocityW[i,j,k+1], [mesh.dz[k]; mesh.dz[k+1]], (mesh.dx[i] * mesh.dy[j])
        )
        A[1, phi.gIndex[i,j,k+1]] = at
        b += b6
    end

    #Center Coefficent
    if phi.bounds[i,j,k]
        @inbounds ac, b0 = _diffusion_centralDifference_central_(
            i,
            j,
            k,
            velocityU,
            velocityV,
            velocityW,
            phi,
            bounds,
            mesh;
            T = T,
        )
        A[1, phi.gIndex[i,j,k]] = ac + aec - awc + anc - asc + atc - abc
        b[phi.gIndex[i,j,k]] += b0
    else
        A[1, phi.gIndex[i,j,k]] = ac + aec - awc + anc - asc + atc - abc
    end

    return A, b
end
