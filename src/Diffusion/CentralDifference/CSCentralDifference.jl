"""

"""
function _discretize_diffusion_centralDifference_ end

function _discretize_diffusion_centralDifference_(
    i::Signed,
    phi::CSPhi1D,
    bounds::Dict{String,BoundsStructured},
    material::CSMaterial1D,
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
    b1 = 0.0
    ae = 0.0
    b2 = 0.0

    #West Coefficents
    if (i != 1) && (phi.onoff[i-1])
        gamma = gamma_interpolation(
            [mesh.dx[i]; mesh.dx[i-1]], [material.Γ[i]; material.Γ[i-1]];
            interpolation = interpolation
        )
        @inbounds aw, b1 = _diffusion_centralDifference_neighbors_(
            gamma, [mesh.dx[i]; mesh.dx[i-1]], (1.0)
        )
        A[1, phi.gIndex[i-1]] = aw
        b += b1
    end

    #East Coefficents
    if (i != mesh.l1)  && (phi.onoff[i+1])
        gamma = gamma_interpolation(
            [mesh.dx[i]; mesh.dx[i+1]], [material.Γ[i]; material.Γ[i+1]];
            interpolation = interpolation
        )
        @inbounds ae, b2 = _diffusion_centralDifference_neighbors_(
            gamma, [mesh.dx[i]; mesh.dx[i+1]], (1.0)
        )
        A[1, phi.gIndex[i+1]] = ae
        b += b2
    end

    #Center Coefficent
    if phi.bounds[i]
        @inbounds ac, b0 = _diffusion_centralDifference_central_(
            i,
            phi,
            bounds,
            mesh;
            T = T,
        )
        A[1, phi.gIndex[i]] = ac - (aw + ae)
        b += b0
    else
        A[1, phi.gIndex[i]] = -1.0 * (aw + ae)
    end

    return A, b
end

function _discretize_diffusion_centralDifference_(
    i::Signed,
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
    b1 = 0.0
    ae = 0.0
    b2 = 0.0

    #West Coefficents
    if (i != 1) && (phi.onoff[i-1])
        gamma = material.Γ
        @inbounds aw, b1 = _diffusion_centralDifference_neighbors_(
            gamma, [mesh.dx[i]; mesh.dx[i-1]], (1.0)
        )
        A[1, phi.gIndex[i-1]] = aw
        b += b1
    end

    #East Coefficents
    if (i != mesh.l1) && (phi.onoff[i+1])
        gamma = material.Γ
        @inbounds ae, b2 = _diffusion_centralDifference_neighbors_(
            gamma, [mesh.dx[i]; mesh.dx[i+1]], (1.0)
        )
        A[1, phi.gIndex[i+1]] = ae
        b += b2
    end

    #Center Coefficent
    if phi.bounds[i]
        @inbounds ac, b0 = _diffusion_centralDifference_central_(
            i,
            phi,
            bounds,
            mesh;
            T = T,
        )
        A[1, phi.gIndex[i]] = ac - (aw + ae)
        b += b0
    else
        A[1, phi.gIndex[i]] = -1.0 * (aw + ae)
    end

    return A, b
end

function _discretize_diffusion_centralDifference_(
    i::Signed,
    j::Signed,
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
    b1 = 0.0
    ae = 0.0
    b2 = 0.0
    as = 0.0
    b3 = 0.0
    an = 0.0
    b4 = 0.0

    #West Coefficents
    if (i != 1) && (phi.onoff[i-1,j])
        gamma = gamma_interpolation(
            [mesh.dx[i]; mesh.dx[i-1]], [material.Γ[i,j]; material.Γ[i-1,j]];
            interpolation = interpolation
        )
        @inbounds aw, b1 = _diffusion_centralDifference_neighbors_(
            gamma, [mesh.dx[i]; mesh.dx[i-1]], (mesh.dy[j])
        )
        A[1, phi.gIndex[i-1,j]] = aw
        b += b1
    end

    #East Coefficents
    if (i != mesh.l1) && (phi.onoff[i+1,j])
        gamma = gamma_interpolation(
            [mesh.dx[i]; mesh.dx[i+1]], [material.Γ[i,j]; material.Γ[i+1,j]];
            interpolation = interpolation
        )
        @inbounds ae, b2 = _diffusion_centralDifference_neighbors_(
            gamma, [mesh.dx[i]; mesh.dx[i+1]], (mesh.dy[j])
        )
        A[1, phi.gIndex[i+1,j]] = ae
        b += b2
    end

    #South Coefficents
    if (j != 1) && (phi.onoff[i,j-1])
        gamma = gamma_interpolation(
            [mesh.dy[j]; mesh.dy[j-1]], [material.Γ[i,j]; material.Γ[i,j-1]];
            interpolation = interpolation
        )
        @inbounds as, b3 = _diffusion_centralDifference_neighbors_(
            gamma, [mesh.dy[j]; mesh.dy[j-1]], (mesh.dx[i])
        )
        A[1, phi.gIndex[i,j-1]] = as
        b += b3
    end

    #North Coefficents
    if (j != mesh.m1) && (phi.onoff[i,j+1])
        gamma = gamma_interpolation(
            [mesh.dy[j]; mesh.dy[j+1]], [material.Γ[i,j]; material.Γ[i,j+1]];
            interpolation = interpolation
        )
        @inbounds an, b4 = _diffusion_centralDifference_neighbors_(
            gamma, [mesh.dy[j]; mesh.dy[j+1]], (mesh.dx[i])
        )
        A[1, phi.gIndex[i,j+1]] = an
        b += b4
    end

    #Center Coefficent
    if phi.bounds[i,j]
        @inbounds ac, b0 = _diffusion_centralDifference_central_(
            i,
            j,
            phi,
            bounds,
            mesh;
            T = T,
        )
        A[1, phi.gIndex[i,j]] = ac - (aw + ae + as + an)
        b += b0
    else
        A[1, phi.gIndex[i,j]] = -1.0 * (aw + ae + as + an)
    end

    return A, b
end

function _discretize_diffusion_centralDifference_(
    i::Signed,
    j::Signed,
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
    b1 = 0.0
    ae = 0.0
    b2 = 0.0
    as = 0.0
    b3 = 0.0
    an = 0.0
    b4 = 0.0

    #West Coefficents
    if (i != 1) && (phi.onoff[i-1,j])
        gamma = material.Γ
        @inbounds aw, b1 = _diffusion_centralDifference_neighbors_(
            gamma, [mesh.dx[i]; mesh.dx[i-1]], (mesh.dy[j])
        )
        A[1, phi.gIndex[i-1,j]] = aw
        b += b1
    end

    #East Coefficents
    if (i != mesh.l1) && (phi.onoff[i+1,j])
        gamma = material.Γ
        @inbounds ae, b2 = _diffusion_centralDifference_neighbors_(
            gamma, [mesh.dx[i]; mesh.dx[i+1]], (mesh.dy[j])
        )
        A[1, phi.gIndex[i+1,j]] = ae
        b += b2
    end

    #South Coefficents
    if (j != 1) && (phi.onoff[i,j-1])
        gamma = material.Γ
        @inbounds as, b3 = _diffusion_centralDifference_neighbors_(
            gamma, [mesh.dy[j]; mesh.dy[j-1]], (mesh.dx[i])
        )
        A[1, phi.gIndex[i,j-1]] = as
        b += b3
    end

    #North Coefficents
    if (j != mesh.m1) && (phi.onoff[i,j+1])
        gamma = material.Γ
        @inbounds an, b4 = _diffusion_centralDifference_neighbors_(
            gamma, [mesh.dy[j]; mesh.dy[j+1]], (mesh.dx[i])
        )
        A[1, phi.gIndex[i,j+1]] = an
        b += b4
    end

    #Center Coefficent
    if phi.bounds[i,j]
        @inbounds ac, b0 = _diffusion_centralDifference_central_(
            i,
            j,
            phi,
            bounds,
            mesh;
            T = T,
        )
        A[1, phi.gIndex[i,j]] = ac - (aw + ae + as + an)
        b += b0
    else
        A[1, phi.gIndex[i,j]] = -1.0 * (aw + ae + as + an)
    end

    return A, b
end

function _discretize_diffusion_centralDifference_(
    i::Signed,
    j::Signed,
    k::Signed,
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
    b1 = 0.0
    ae = 0.0
    b2 = 0.0
    as = 0.0
    b3 = 0.0
    an = 0.0
    b4 = 0.0
    ab = 0.0
    b5 = 0.0
    at = 0.0
    b6 = 0.0

    #West Coefficents
    if (i != 1) && (phi.onoff[i-1,j,k])
        gamma = gamma_interpolation(
            [mesh.dx[i]; mesh.dx[i-1]], [material.Γ[i,j,k]; material.Γ[i-1,j,k]];
            interpolation = interpolation
        )
        @inbounds aw, b1 = _diffusion_centralDifference_neighbors_(
            gamma, [mesh.dx[i]; mesh.dx[i-1]], (mesh.dy[j] * mesh.dz[k])
        )
        A[1, phi.gIndex[i-1,j,k]] = aw
        b += b1
    end

    #East Coefficents
    if (i != mesh.l1) && (phi.onoff[i+1,j,k])
        gamma = gamma_interpolation(
            [mesh.dx[i]; mesh.dx[i+1]], [material.Γ[i,j,k]; material.Γ[i+1,j,k]];
            interpolation = interpolation
        )
        @inbounds ae, b2 = _diffusion_centralDifference_neighbors_(
            gamma, [mesh.dx[i]; mesh.dx[i+1]], (mesh.dy[j] * mesh.dz[k])
        )
        A[1, phi.gIndex[i+1,j,k]] = ae
        b += b2
    end

    #South Coefficents
    if (j != 1) && (phi.onoff[i,j-1,k])
        gamma = gamma_interpolation(
            [mesh.dy[j]; mesh.dy[j-1]], [material.Γ[i,j,k]; material.Γ[i,j-1,k]];
            interpolation = interpolation
        )
        @inbounds as, b3 = _diffusion_centralDifference_neighbors_(
            gamma, [mesh.dy[j]; mesh.dy[j-1]], (mesh.dx[i] * mesh.dz[k])
        )
        A[1, phi.gIndex[i,j-1,k]] = as
        b += b3
    end

    #North Coefficents
    if (j != mesh.m1) && (phi.onoff[i,j+1,k])
        gamma = gamma_interpolation(
            [mesh.dy[j]; mesh.dy[j+1]], [material.Γ[i,j,k]; material.Γ[i,j+1,k]];
            interpolation = interpolation
        )
        @inbounds an, b4 = _diffusion_centralDifference_neighbors_(
            gamma, [mesh.dy[j]; mesh.dy[j+1]], (mesh.dx[i] * mesh.dz[k])
        )
        A[1, phi.gIndex[i,j+1,k]] = an
        b += b4
    end

    #Bottom Coefficents
    if (k != 1) && (phi.onoff[i,j,k-1])
        gamma = gamma_interpolation(
            [mesh.dz[k]; mesh.dz[k-1]], [material.Γ[i,j,k]; material.Γ[i,j,k-1]];
            interpolation = interpolation
        )
        @inbounds ab, b5 = _diffusion_centralDifference_neighbors_(
            gamma, [mesh.dz[k]; mesh.dz[k-1]], (mesh.dx[i] * mesh.dy[j])
        )
        A[1, phi.gIndex[i,j,k-1]] = ab
        b += b5
    end

    #Top Coefficents
    if (k != mesh.n1) && (phi.onoff[i,j,k+1])
        gamma = gamma_interpolation(
            [mesh.dz[k]; mesh.dz[k+1]], [material.Γ[i,j,k]; material.Γ[i,j,k+1]];
            interpolation = interpolation
        )
        @inbounds at, b6 = _diffusion_centralDifference_neighbors_(
            gamma, [mesh.dz[k]; mesh.dz[k+1]], (mesh.dx[i] * mesh.dy[j])
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
            phi,
            bounds,
            mesh;
            T = T,
        )
        A[1, phi.gIndex[i,j,k]] = ac - (aw + ae + as + an + ab + at)
        b += b0
    else
        A[1, phi.gIndex[i,j,k]] = -1.0 * (aw + ae + as + an + ab + at)
    end

    return A, b
end

function _discretize_diffusion_centralDifference_(
    i::Signed,
    j::Signed,
    k::Signed,
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
    b1 = 0.0
    ae = 0.0
    b2 = 0.0
    as = 0.0
    b3 = 0.0
    an = 0.0
    b4 = 0.0
    ab = 0.0
    b5 = 0.0
    at = 0.0
    b6 = 0.0

    #West Coefficents
    if (i != 1) && (phi.onoff[i-1,j,k])
        gamma = material.Γ
        @inbounds aw, b1 = _diffusion_centralDifference_neighbors_(
            gamma, [mesh.dx[i]; mesh.dx[i-1]], (mesh.dy[j] * mesh.dz[k])
        )
        A[1, phi.gIndex[i-1,j,k]] = aw
        b += b1
    end

    #East Coefficents
    if (i != mesh.l1) && (phi.onoff[i+1,j,k])
        gamma = material.Γ
        @inbounds ae, b2 = _diffusion_centralDifference_neighbors_(
            gamma, [mesh.dx[i]; mesh.dx[i+1]], (mesh.dy[j] * mesh.dz[k])
        )
        A[1, phi.gIndex[i+1,j,k]] = ae
        b += b2
    end

    #South Coefficents
    if (j != 1) && (phi.onoff[i,j-1,k])
        gamma = material.Γ
        @inbounds as, b3 = _diffusion_centralDifference_neighbors_(
            gamma, [mesh.dy[j]; mesh.dy[j-1]], (mesh.dx[i] * mesh.dz[k])
        )
        A[1, phi.gIndex[i,j-1,k]] = as
        b += b3
    end

    #North Coefficents
    if (j != mesh.m1) && (phi.onoff[i,j+1,k])
        gamma = material.Γ
        @inbounds an, b4 = _diffusion_centralDifference_neighbors_(
            gamma, [mesh.dy[j]; mesh.dy[j+1]], (mesh.dx[i] * mesh.dz[k])
        )
        A[1, phi.gIndex[i,j+1,k]] = an
        b += b4
    end

    #Bottom Coefficents
    if (k != 1) && (phi.onoff[i,j,k-1])
        gamma = material.Γ
        @inbounds ab, b5 = _diffusion_centralDifference_neighbors_(
            gamma, [mesh.dz[k]; mesh.dz[k-1]], (mesh.dx[i] * mesh.dy[j])
        )
        A[1, phi.gIndex[i,j,k-1]] = ab
        b += b5
    end

    #Top Coefficents
    if (k != mesh.n1) && (phi.onoff[i,j,k+1])
        gamma = material.Γ
        @inbounds at, b6 = _diffusion_centralDifference_neighbors_(
            gamma, [mesh.dz[k]; mesh.dz[k+1]], (mesh.dx[i] * mesh.dy[j])
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
            phi,
            bounds,
            mesh;
            T = T,
        )
        A[1, phi.gIndex[i,j,k]] = ac - (aw + ae + as + an + ab + at)
        b += b0
    else
        A[1, phi.gIndex[i,j,k]] = -1.0 * (aw + ae + as + an + ab + at)
    end

    return A, b
end
