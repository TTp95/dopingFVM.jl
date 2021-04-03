"""

"""
function _discretize_diffusion_centralDifference_ end

function _discretize_diffusion_centralDifference_(
    phi::CSPhi1D,
    bounds::Dict{String,BoundsStructured},
    material::CSMaterial1D,
    mesh::UnionCSMesh1D;
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
    mthreads::Bool = false,
    sparrays::Bool = true,
    interpolation::Signed = 2,
)
    n_equations = maximum_globalIndex(phi)

    if mthreads

    elseif sparrays
        n = 0
        AI = zeros(N, (3 * n_equations))
        AJ = zeros(N, (3 * n_equations))
        AV = zeros(T, (3 * n_equations))
    elseif !sparrays
        A = zeros(T, n_equations, n_equations)
    end

    b = zeros(T, n_equations)

    if mthreads
        #Base.Threads.@threads
    elseif sparrays
        for i in 1:mesh.l1
            if phi.onoff[i]
                id = phi.gIndex[i]

                #Auxiliar variables
                ac = 0.0
                b0 = 0.0
                aw = 0.0
                b1 = 0.0
                ae = 0.0
                b2 = 0.0

                #West Coefficents
                if (i != 1) && (mesh.l1 != 1) && (phi.onoff[i-1])
                    gamma = gamma_interpolation(
                        mesh.dx[i], mesh.dx[i-1], material.Γ[i], material.Γ[i-1];
                        interpolation = interpolation
                    )
                    @inbounds aw, b1 = _diffusion_centralDifference_neighbors_(
                        gamma, mesh.dx[i], mesh.dx[i-1], (1.0)
                    )
                    n += 1
                    AI[n] = id
                    AJ[n] = phi.gIndex[i-1]
                    AV[n] = aw
                    b[id] += b1
                end

                #East Coefficents
                if (i != mesh.l1) && (mesh.l1 != 1)  && (phi.onoff[i+1])
                    gamma = gamma_interpolation(
                        mesh.dx[i], mesh.dx[i+1], material.Γ[i], material.Γ[i+1];
                        interpolation = interpolation
                    )
                    @inbounds ae, b2 = _diffusion_centralDifference_neighbors_(
                        gamma, mesh.dx[i], mesh.dx[i+1], (1.0)
                    )
                    n += 1
                    AI[n] = id
                    AJ[n] = phi.gIndex[i+1]
                    AV[n] = ae
                    b[id] += b2
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
                    n += 1
                    AI[n] = id
                    AJ[n] = id
                    AV[n] = ac - (aw + ae)
                    b[id] += b0
                else
                    n += 1
                    AI[n] = id
                    AJ[n] = id
                    AV[n] = -1.0 * (aw + ae)
                end
            end
        end
    elseif !sparrays
        for i in 1:mesh.l1
            if phi.onoff[i]
                id = phi.gIndex[i]

                #Auxiliar variables
                ac = 0.0
                b0 = 0.0
                aw = 0.0
                b1 = 0.0
                ae = 0.0
                b2 = 0.0

                #West Coefficents
                if (i != 1) && (mesh.l1 != 1) && (phi.onoff[i-1])
                    gamma = gamma_interpolation(
                        mesh.dx[i], mesh.dx[i-1], material.Γ[i], material.Γ[i-1];
                        interpolation = interpolation
                    )
                    @inbounds aw, b1 = _diffusion_centralDifference_neighbors_(
                        gamma, mesh.dx[i], mesh.dx[i-1], (1.0)
                    )
                    A[id, phi.gIndex[i-1]] = aw
                    b[id] += b1
                end

                #East Coefficents
                if (i != mesh.l1) && (mesh.l1 != 1)  && (phi.onoff[i+1])
                    gamma = gamma_interpolation(
                        mesh.dx[i], mesh.dx[i+1], material.Γ[i], material.Γ[i+1];
                        interpolation = interpolation
                    )
                    @inbounds ae, b2 = _diffusion_centralDifference_neighbors_(
                        gamma, mesh.dx[i], mesh.dx[i+1], (1.0)
                    )
                    A[id, phi.gIndex[i+1]] = ae
                    b[id] += b2
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
                    A[id, phi.gIndex[i]] = ac - (aw + ae)
                    b[id] += b0
                else
                    A[id, phi.gIndex[i]] = -1.0 * (aw + ae)
                end
            end
        end
    end

    if sparrays
        A = sparse(AI[1:n], AJ[1:n], AV[1:n])
    end

    return A, b
end

function _discretize_diffusion_centralDifference_(
    phi::CSPhi1D,
    bounds::Dict{String,BoundsStructured},
    material::UnionCSConstantMaterial,
    mesh::UnionCSMesh1D;
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
    mthreads::Bool = false,
    sparrays::Bool = true,
    interpolation::Signed = 2,
)
    n_equations = maximum_globalIndex(phi)

    if mthreads

    elseif sparrays
        n = 0
        AI = zeros(N, (3 * n_equations))
        AJ = zeros(N, (3 * n_equations))
        AV = zeros(T, (3 * n_equations))
    elseif !sparrays
        A = zeros(T, n_equations, n_equations)
    end

    b = zeros(T, n_equations)

    if mthreads
        #Base.Threads.@threads
    elseif sparrays
        for i in 1:mesh.l1
            if phi.onoff[i]
                id = phi.gIndex[i]

                #Auxiliar variables
                ac = 0.0
                b0 = 0.0
                aw = 0.0
                b1 = 0.0
                ae = 0.0
                b2 = 0.0

                #West Coefficents
                if (i != 1) && (mesh.l1 != 1) && (phi.onoff[i-1])
                    gamma = material.Γ
                    @inbounds aw, b1 = _diffusion_centralDifference_neighbors_(
                        gamma, mesh.dx[i], mesh.dx[i-1], (1.0)
                    )
                    n += 1
                    AI[n] = id
                    AJ[n] = phi.gIndex[i-1]
                    AV[n] = aw
                    b[id] += b1
                end

                #East Coefficents
                if (i != mesh.l1) && (mesh.l1 != 1)  && (phi.onoff[i+1])
                    gamma = material.Γ
                    @inbounds ae, b2 = _diffusion_centralDifference_neighbors_(
                        gamma, mesh.dx[i], mesh.dx[i+1], (1.0)
                    )
                    n += 1
                    AI[n] = id
                    AJ[n] = phi.gIndex[i+1]
                    AV[n] = ae
                    b[id] += b2
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
                    n += 1
                    AI[n] = id
                    AJ[n] = id
                    AV[n] = ac - (aw + ae)
                    b[id] += b0
                else
                    n += 1
                    AI[n] = id
                    AJ[n] = id
                    AV[n] = -1.0 * (aw + ae)
                end
            end
        end
    elseif !sparrays
        for i in 1:mesh.l1
            if phi.onoff[i]
                id = phi.gIndex[i]

                #Auxiliar variables
                ac = 0.0
                b0 = 0.0
                aw = 0.0
                b1 = 0.0
                ae = 0.0
                b2 = 0.0

                #West Coefficents
                if (i != 1) && (mesh.l1 != 1) && (phi.onoff[i-1])
                    gamma = material.Γ
                    @inbounds aw, b1 = _diffusion_centralDifference_neighbors_(
                        gamma, mesh.dx[i], mesh.dx[i-1], (1.0)
                    )
                    A[id, phi.gIndex[i-1]] = aw
                    b[id] += b1
                end

                #East Coefficents
                if (i != mesh.l1) && (mesh.l1 != 1)  && (phi.onoff[i+1])
                    gamma = material.Γ
                    @inbounds ae, b2 = _diffusion_centralDifference_neighbors_(
                        gamma, mesh.dx[i], mesh.dx[i+1], (1.0)
                    )
                    A[id, phi.gIndex[i+1]] = ae
                    b[id] += b2
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
                    A[id, phi.gIndex[i]] = ac - (aw + ae)
                    b[id] += b0
                else
                    A[id, phi.gIndex[i]] = -1.0 * (aw + ae)
                end
            end
        end
    end

    if sparrays
        A = sparse(AI[1:n], AJ[1:n], AV[1:n])
    end

    return A, b
end

function _discretize_diffusion_centralDifference_(
    phi::CSPhi2D,
    bounds::Dict{String,BoundsStructured},
    material::CSMaterial2D,
    mesh::UnionCSMesh2D;
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
    mthreads::Bool = false,
    sparrays::Bool = true,
    interpolation::Signed = 2,
)
    n_equations = maximum_globalIndex(phi)

    if mthreads

    elseif sparrays
        n = 0
        AI = zeros(N, (5 * n_equations))
        AJ = zeros(N, (5 * n_equations))
        AV = zeros(T, (5 * n_equations))
    elseif !sparrays
        A = zeros(T, n_equations, n_equations)
    end

    b = zeros(T, n_equations)

    if mthreads
        # Base.Threads.@threads
    elseif sparrays
        for i in 1:mesh.l1
            for j in 1:mesh.m1
                if phi.onoff[i,j]
                    id = phi.gIndex[i,j]

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
                    if (i != 1) && (mesh.l1 != 1) && (phi.onoff[i-1,j])
                        gamma = gamma_interpolation(
                            mesh.dx[i], mesh.dx[i-1], material.Γ[i,j], material.Γ[i-1,j];
                            interpolation = interpolation
                        )
                        @inbounds aw, b1 = _diffusion_centralDifference_neighbors_(
                            gamma, mesh.dx[i], mesh.dx[i-1], (mesh.dy[j])
                        )
                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i-1,j]
                        AV[n] = aw
                        b[id] += b1
                    end

                    #East Coefficents
                    if (i != mesh.l1) && (mesh.l1 != 1) && (phi.onoff[i+1,j])
                        gamma = gamma_interpolation(
                            mesh.dx[i], mesh.dx[i+1], material.Γ[i,j], material.Γ[i+1,j];
                            interpolation = interpolation
                        )
                        @inbounds ae, b2 = _diffusion_centralDifference_neighbors_(
                            gamma, mesh.dx[i], mesh.dx[i+1], (mesh.dy[j])
                        )
                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i+1,j]
                        AV[n] = ae
                        b[id] += b2
                    end

                    #South Coefficents
                    if (j != 1) && (mesh.m1 != 1) && (phi.onoff[i,j-1])
                        gamma = gamma_interpolation(
                            mesh.dy[j], mesh.dy[j-1], material.Γ[i,j], material.Γ[i,j-1];
                            interpolation = interpolation
                        )
                        @inbounds as, b3 = _diffusion_centralDifference_neighbors_(
                            gamma, mesh.dy[j], mesh.dy[j-1], (mesh.dx[i])
                        )
                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i,j-1]
                        AV[n] = as
                        b[id] += b3
                    end

                    #North Coefficents
                    if (j != mesh.m1) && (mesh.m1 != 1) && (phi.onoff[i,j+1])
                        gamma = gamma_interpolation(
                            mesh.dy[j], mesh.dy[j+1], material.Γ[i,j], material.Γ[i,j+1];
                            interpolation = interpolation
                        )
                        @inbounds an, b4 = _diffusion_centralDifference_neighbors_(
                            gamma, mesh.dy[j], mesh.dy[j+1], (mesh.dx[i])
                        )
                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i,j+1]
                        AV[n] = an
                        b[id] += b4
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
                        n += 1
                        AI[n] = id
                        AJ[n] = id
                        AV[n] = ac - (aw + ae + as + an)
                        b[id] += b0
                    else
                        n += 1
                        AI[n] = id
                        AJ[n] = id
                        AV[n] = -1.0 * (aw + ae + as + an)
                    end
                end
            end
        end
    elseif !sparrays
        for i in 1:mesh.l1
            for j in 1:mesh.m1
                if phi.onoff[i,j]
                    id = phi.gIndex[i,j]

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
                    if (i != 1) && (mesh.l1 != 1) && (phi.onoff[i-1,j])
                        gamma = gamma_interpolation(
                            mesh.dx[i], mesh.dx[i-1], material.Γ[i,j], material.Γ[i-1,j];
                            interpolation = interpolation
                        )
                        @inbounds aw, b1 = _diffusion_centralDifference_neighbors_(
                            gamma, mesh.dx[i], mesh.dx[i-1], (mesh.dy[j])
                        )
                        A[id, phi.gIndex[i-1,j]] = aw
                        b[id] += b1
                    end

                    #East Coefficents
                    if (i != mesh.l1) && (mesh.l1 != 1) && (phi.onoff[i+1,j])
                        gamma = gamma_interpolation(
                            mesh.dx[i], mesh.dx[i+1], material.Γ[i,j], material.Γ[i+1,j];
                            interpolation = interpolation
                        )
                        @inbounds ae, b2 = _diffusion_centralDifference_neighbors_(
                            gamma, mesh.dx[i], mesh.dx[i+1], (mesh.dy[j])
                        )
                        A[id, phi.gIndex[i+1,j]] = ae
                        b[id] += b2
                    end

                    #South Coefficents
                    if (j != 1) && (mesh.m1 != 1) && (phi.onoff[i,j-1])
                        gamma = gamma_interpolation(
                            mesh.dy[j], mesh.dy[j-1], material.Γ[i,j], material.Γ[i,j-1];
                            interpolation = interpolation
                        )
                        @inbounds as, b3 = _diffusion_centralDifference_neighbors_(
                            gamma, mesh.dy[j], mesh.dy[j-1], (mesh.dx[i])
                        )
                        A[id, phi.gIndex[i,j-1]] = as
                        b[id] += b3
                    end

                    #North Coefficents
                    if (j != mesh.m1) && (mesh.m1 != 1) && (phi.onoff[i,j+1])
                        gamma = gamma_interpolation(
                            mesh.dy[j], mesh.dy[j+1], material.Γ[i,j], material.Γ[i,j+1];
                            interpolation = interpolation
                        )
                        @inbounds an, b4 = _diffusion_centralDifference_neighbors_(
                            gamma, mesh.dy[j], mesh.dy[j+1], (mesh.dx[i])
                        )
                        A[id, phi.gIndex[i,j+1]] = an
                        b[id] += b4
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
                        A[id, id] = ac - (aw + ae + as + an)
                        b[id] += b0
                    else
                        A[id, id] = -1.0 * (aw + ae + as + an)
                    end

                end
            end
        end
    end

    if sparrays
        A = sparse(AI[1:n], AJ[1:n], AV[1:n])
    end

    return A, b
end

function _discretize_diffusion_centralDifference_(
    phi::CSPhi2D,
    bounds::Dict{String,BoundsStructured},
    material::UnionCSConstantMaterial,
    mesh::UnionCSMesh2D;
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
    mthreads::Bool = false,
    sparrays::Bool = true,
    interpolation::Signed = 2,
)
    n_equations = maximum_globalIndex(phi)

    if mthreads

    elseif sparrays
        n = 0
        AI = zeros(N, (5 * n_equations))
        AJ = zeros(N, (5 * n_equations))
        AV = zeros(T, (5 * n_equations))
    elseif !sparrays
        A = zeros(T, n_equations, n_equations)
    end

    b = zeros(T, n_equations)

    if mthreads
        # Base.Threads.@threads
    elseif sparrays
        for i in 1:mesh.l1
            for j in 1:mesh.m1
                if phi.onoff[i,j]
                    id = phi.gIndex[i,j]

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
                    if (i != 1) && (mesh.l1 != 1) && (phi.onoff[i-1,j])
                        gamma = material.Γ
                        @inbounds aw, b1 = _diffusion_centralDifference_neighbors_(
                            gamma, mesh.dx[i], mesh.dx[i-1], (mesh.dy[j])
                        )
                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i-1,j]
                        AV[n] = aw
                        b[id] += b1
                    end

                    #East Coefficents
                    if (i != mesh.l1) && (mesh.l1 != 1) && (phi.onoff[i+1,j])
                        gamma = material.Γ
                        @inbounds ae, b2 = _diffusion_centralDifference_neighbors_(
                            gamma, mesh.dx[i], mesh.dx[i+1], (mesh.dy[j])
                        )
                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i+1,j]
                        AV[n] = ae
                        b[id] += b2
                    end

                    #South Coefficents
                    if (j != 1) && (mesh.m1 != 1) && (phi.onoff[i,j-1])
                        gamma = material.Γ
                        @inbounds as, b3 = _diffusion_centralDifference_neighbors_(
                            gamma, mesh.dy[j], mesh.dy[j-1], (mesh.dx[i])
                        )
                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i,j-1]
                        AV[n] = as
                        b[id] += b3
                    end

                    #North Coefficents
                    if (j != mesh.m1) && (mesh.m1 != 1) && (phi.onoff[i,j+1])
                        gamma = material.Γ
                        @inbounds an, b4 = _diffusion_centralDifference_neighbors_(
                            gamma, mesh.dy[j], mesh.dy[j+1], (mesh.dx[i])
                        )
                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i,j+1]
                        AV[n] = an
                        b[id] += b4
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
                        n += 1
                        AI[n] = id
                        AJ[n] = id
                        AV[n] = ac - (aw + ae + as + an)
                        b[id] += b0
                    else
                        n += 1
                        AI[n] = id
                        AJ[n] = id
                        AV[n] = -1.0 * (aw + ae + as + an)
                    end

                end
            end
        end
    elseif !sparrays
        for i in 1:mesh.l1
            for j in 1:mesh.m1
                if phi.onoff[i,j]
                    id = phi.gIndex[i,j]

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
                    if (i != 1) && (mesh.l1 != 1) && (phi.onoff[i-1,j])
                        gamma = material.Γ
                        @inbounds aw, b1 = _diffusion_centralDifference_neighbors_(
                            gamma, mesh.dx[i], mesh.dx[i-1], (mesh.dy[j])
                        )
                        A[id, phi.gIndex[i-1,j]] = aw
                        b[id] += b1
                    end

                    #East Coefficents
                    if (i != mesh.l1) && (mesh.l1 != 1) && (phi.onoff[i+1,j])
                        gamma = material.Γ
                        @inbounds ae, b2 = _diffusion_centralDifference_neighbors_(
                            gamma, mesh.dx[i], mesh.dx[i+1], (mesh.dy[j])
                        )
                        A[id, phi.gIndex[i+1,j]] = ae
                        b[id] += b2
                    end

                    #South Coefficents
                    if (j != 1) && (mesh.m1 != 1) && (phi.onoff[i,j-1])
                        gamma = material.Γ
                        @inbounds as, b3 = _diffusion_centralDifference_neighbors_(
                            gamma, mesh.dy[j], mesh.dy[j-1], (mesh.dx[i])
                        )
                        A[id, phi.gIndex[i,j-1]] = as
                        b[id] += b3
                    end

                    #North Coefficents
                    if (j != mesh.m1) && (mesh.m1 != 1) && (phi.onoff[i,j+1])
                        gamma = material.Γ
                        @inbounds an, b4 = _diffusion_centralDifference_neighbors_(
                            gamma, mesh.dy[j], mesh.dy[j+1], (mesh.dx[i])
                        )
                        A[id, phi.gIndex[i,j+1]] = an
                        b[id] += b4
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
                        A[id, id] = ac - (aw + ae + as + an)
                        b[id] += b0
                    else
                        A[id, id] = -1.0 * (aw + ae + as + an)
                    end

                end
            end
        end
    end

    if sparrays
        A = sparse(AI[1:n], AJ[1:n], AV[1:n])
    end

    return A, b
end

function _discretize_diffusion_centralDifference_(
    phi::CSPhi3D,
    bounds::Dict{String,BoundsStructured},
    material::CSMaterial3D,
    mesh::UnionCSMesh3D;
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
    mthreads::Bool = false,
    sparrays::Bool = true,
    interpolation::Signed = 2,
)
    n_equations = maximum_globalIndex(phi)

    if mthreads

    elseif sparrays
        n = 0
        AI = zeros(N, (7 * n_equations))
        AJ = zeros(N, (7 * n_equations))
        AV = zeros(T, (7 * n_equations))
    elseif !sparrays
        A = zeros(T, n_equations, n_equations)
    end

    b = zeros(T, n_equations)

    if mthreads
        # Base.Threads.@threads
    elseif sparrays
        for i in 1:mesh.l1
            for j in 1:mesh.m1
                for k in 1:mesh.n1
                    if phi.onoff[i,j,k]
                        id = phi.gIndex[i,j,k]

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
                        if (i != 1) && (mesh.l1 != 1) && (phi.onoff[i-1,j,k])
                            gamma = gamma_interpolation(
                                mesh.dx[i], mesh.dx[i-1], material.Γ[i,j,k], material.Γ[i-1,j,k];
                                interpolation = interpolation
                            )
                            @inbounds aw, b1 = _diffusion_centralDifference_neighbors_(
                                gamma, mesh.dx[i], mesh.dx[i-1], (mesh.dy[j] * mesh.dz[k])
                            )
                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i-1,j,k]
                            AV[n] = aw
                            b[id] += b1
                        end

                        #East Coefficents
                        if (i != mesh.l1) && (mesh.l1 != 1) && (phi.onoff[i+1,j,k])
                            gamma = gamma_interpolation(
                                mesh.dx[i], mesh.dx[i+1], material.Γ[i,j,k], material.Γ[i+1,j,k];
                                interpolation = interpolation
                            )
                            @inbounds ae, b2 = _diffusion_centralDifference_neighbors_(
                                gamma, mesh.dx[i], mesh.dx[i+1], (mesh.dy[j] * mesh.dz[k])
                            )
                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i+1,j,k]
                            AV[n] = ae
                            b[id] += b2
                        end

                        #South Coefficents
                        if (j != 1) && (mesh.m1 != 1) && (phi.onoff[i,j-1,k])
                            gamma = gamma_interpolation(
                                mesh.dy[j], mesh.dy[j-1], material.Γ[i,j,k], material.Γ[i,j-1,k];
                                interpolation = interpolation
                            )
                            @inbounds as, b3 = _diffusion_centralDifference_neighbors_(
                                gamma, mesh.dy[j], mesh.dy[j-1], (mesh.dx[i] * mesh.dz[k])
                            )
                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i,j-1,k]
                            AV[n] = as
                            b[id] += b3
                        end

                        #North Coefficents
                        if (j != mesh.m1) && (mesh.m1 != 1) && (phi.onoff[i,j+1,k])
                            gamma = gamma_interpolation(
                                mesh.dy[j], mesh.dy[j+1], material.Γ[i,j,k], material.Γ[i,j+1,k];
                                interpolation = interpolation
                            )
                            @inbounds an, b4 = _diffusion_centralDifference_neighbors_(
                                gamma, mesh.dy[j], mesh.dy[j+1], (mesh.dx[i] * mesh.dz[k])
                            )
                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i,j+1,k]
                            AV[n] = an
                            b[id] += b4
                        end

                        #Bottom Coefficents
                        if (k != 1) && (mesh.n1 != 1) && (phi.onoff[i,j,k-1])
                            gamma = gamma_interpolation(
                                mesh.dz[k], mesh.dz[k-1], material.Γ[i,j,k], material.Γ[i,j,k-1];
                                interpolation = interpolation
                            )
                            @inbounds ab, b5 = _diffusion_centralDifference_neighbors_(
                                gamma, mesh.dz[k], mesh.dz[k-1], (mesh.dx[i] * mesh.dy[j])
                            )
                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i,j,k-1]
                            AV[n] = ab
                            b[id] += b5
                        end

                        #Top Coefficents
                        if (k != mesh.n1) && (mesh.n1 != 1) && (phi.onoff[i,j,k+1])
                            gamma = gamma_interpolation(
                                mesh.dz[k], mesh.dz[k+1], material.Γ[i,j,k], material.Γ[i,j,k+1];
                                interpolation = interpolation
                            )
                            @inbounds at, b6 = _diffusion_centralDifference_neighbors_(
                                gamma, mesh.dz[k], mesh.dz[k+1], (mesh.dx[i] * mesh.dy[j])
                            )
                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i,j,k+1]
                            AV[n] = at
                            push!(AI, id)
                            push!(AJ, phi.gIndex[i,j,k+1])
                            push!(AV, at)
                            b[id] += b6
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
                            n += 1
                            AI[n] = id
                            AJ[n] = id
                            AV[n] = ac - (aw + ae + as + an + ab + at)
                            b[id] += b0
                        else
                            n += 1
                            AI[n] = id
                            AJ[n] = id
                            AV[n] = -1.0 * (aw + ae + as + an + ab + at)
                        end
                    end
                end
            end
        end
    elseif !sparrays
        for i in 1:mesh.l1
            for j in 1:mesh.m1
                for k in 1:mesh.n1
                    if phi.onoff[i,j,k]
                        id = phi.gIndex[i,j,k]

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
                        if (i != 1) && (mesh.l1 != 1) && (phi.onoff[i-1,j,k])
                            gamma = gamma_interpolation(
                                mesh.dx[i], mesh.dx[i-1], material.Γ[i,j,k], material.Γ[i-1,j,k];
                                interpolation = interpolation
                            )
                            @inbounds aw, b1 = _diffusion_centralDifference_neighbors_(
                                gamma, mesh.dx[i], mesh.dx[i-1], (mesh.dy[j] * mesh.dz[k])
                            )
                            A[id, phi.gIndex[i-1,j,k]] = aw
                            b[id] += b1
                        end

                        #East Coefficents
                        if (i != mesh.l1) && (mesh.l1 != 1) && (phi.onoff[i+1,j,k])
                            gamma = gamma_interpolation(
                                mesh.dx[i], mesh.dx[i+1], material.Γ[i,j,k], material.Γ[i+1,j,k];
                                interpolation = interpolation
                            )
                            @inbounds ae, b2 = _diffusion_centralDifference_neighbors_(
                                gamma, mesh.dx[i], mesh.dx[i+1], (mesh.dy[j] * mesh.dz[k])
                            )
                            A[id, phi.gIndex[i+1,j,k]] = ae
                            b[id] += b2
                        end

                        #South Coefficents
                        if (j != 1) && (mesh.m1 != 1) && (phi.onoff[i,j-1,k])
                            gamma = gamma_interpolation(
                                mesh.dy[j], mesh.dy[j-1], material.Γ[i,j,k], material.Γ[i,j-1,k];
                                interpolation = interpolation
                            )
                            @inbounds as, b3 = _diffusion_centralDifference_neighbors_(
                                gamma, mesh.dy[j], mesh.dy[j-1], (mesh.dx[i] * mesh.dz[k])
                            )
                            A[id, phi.gIndex[i,j-1,k]] = as
                            b[id] += b3
                        end

                        #North Coefficents
                        if (j != mesh.m1) && (mesh.m1 != 1) && (phi.onoff[i,j+1,k])
                            gamma = gamma_interpolation(
                                mesh.dy[j], mesh.dy[j+1], material.Γ[i,j,k], material.Γ[i,j+1,k];
                                interpolation = interpolation
                            )
                            @inbounds an, b4 = _diffusion_centralDifference_neighbors_(
                                gamma, mesh.dy[j], mesh.dy[j+1], (mesh.dx[i] * mesh.dz[k])
                            )
                            A[id, phi.gIndex[i,j+1,k]] = an
                            b[id] += b4
                        end

                        #Bottom Coefficents
                        if (k != 1) && (mesh.n1 != 1) && (phi.onoff[i,j,k-1])
                            gamma = gamma_interpolation(
                                mesh.dz[k], mesh.dz[k-1], material.Γ[i,j,k], material.Γ[i,j,k-1];
                                interpolation = interpolation
                            )
                            @inbounds ab, b5 = _diffusion_centralDifference_neighbors_(
                                gamma, mesh.dz[k], mesh.dz[k-1], (mesh.dx[i] * mesh.dy[j])
                            )
                            A[id, phi.gIndex[i,j,k-1]] = ab
                            b[id] += b5
                        end

                        #Top Coefficents
                        if (k != mesh.n1) && (mesh.n1 != 1) && (phi.onoff[i,j,k+1])
                            gamma = gamma_interpolation(
                                mesh.dz[k], mesh.dz[k+1], material.Γ[i,j,k], material.Γ[i,j,k+1];
                                interpolation = interpolation
                            )
                            @inbounds at, b6 = _diffusion_centralDifference_neighbors_(
                                gamma, mesh.dz[k], mesh.dz[k+1], (mesh.dx[i] * mesh.dy[j])
                            )
                            A[id, phi.gIndex[i,j,k+1]] = at
                            b[id] += b6
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
                            A[id, id] = ac - (aw + ae + as + an + ab + at)
                            b[id] += b0
                        else
                            A[id, id] = -1.0 * (aw + ae + as + an + ab + at)
                        end
                    end
                end
            end
        end
    end

    if sparrays
        A = sparse(AI[1:n], AJ[1:n], AV[1:n])
    end

    return A, b
end

function _discretize_diffusion_centralDifference_(
    phi::CSPhi3D,
    bounds::Dict{String,BoundsStructured},
    material::UnionCSConstantMaterial,
    mesh::UnionCSMesh3D;
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
    mthreads::Bool = false,
    sparrays::Bool = true,
    interpolation::Signed = 2,
)
    n_equations = maximum_globalIndex(phi)

    if mthreads

    elseif sparrays
        n = 0
        AI = zeros(N, (7 * n_equations))
        AJ = zeros(N, (7 * n_equations))
        AV = zeros(T, (7 * n_equations))
    elseif !sparrays
        A = zeros(T, n_equations, n_equations)
    end

    b = zeros(T, n_equations)

    if mthreads
        # Base.Threads.@threads
    elseif sparrays
        for i in 1:mesh.l1
            for j in 1:mesh.m1
                for k in 1:mesh.n1
                    if phi.onoff[i,j,k]
                        id = phi.gIndex[i,j,k]

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
                        if (i != 1) && (mesh.l1 != 1) && (phi.onoff[i-1,j,k])
                            gamma = material.Γ
                            @inbounds aw, b1 = _diffusion_centralDifference_neighbors_(
                                gamma, mesh.dx[i], mesh.dx[i-1], (mesh.dy[j] * mesh.dz[k])
                            )
                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i,j,k]
                            AV[n] = aw
                            b[id] += b1
                        end

                        #East Coefficents
                        if (i != mesh.l1) && (mesh.l1 != 1) && (phi.onoff[i+1,j,k])
                            gamma = material.Γ
                            @inbounds ae, b2 = _diffusion_centralDifference_neighbors_(
                                gamma, mesh.dx[i], mesh.dx[i+1], (mesh.dy[j] * mesh.dz[k])
                            )
                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i+1,j,k]
                            AV[n] = ae
                            b[id] += b2
                        end

                        #South Coefficents
                        if (j != 1) && (mesh.m1 != 1) && (phi.onoff[i,j-1,k])
                            gamma = material.Γ
                            @inbounds as, b3 = _diffusion_centralDifference_neighbors_(
                                gamma, mesh.dy[j], mesh.dy[j-1], (mesh.dx[i] * mesh.dz[k])
                            )
                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i,j-1,k]
                            AV[n] = as
                            b[id] += b3
                        end

                        #North Coefficents
                        if (j != mesh.m1) && (mesh.m1 != 1) && (phi.onoff[i,j+1,k])
                            gamma = material.Γ
                            @inbounds an, b4 = _diffusion_centralDifference_neighbors_(
                                gamma, mesh.dy[j], mesh.dy[j+1], (mesh.dx[i] * mesh.dz[k])
                            )
                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i,j+1,k]
                            AV[n] = an
                            b[id] += b4
                        end

                        #Bottom Coefficents
                        if (k != 1) && (mesh.n1 != 1) && (phi.onoff[i,j,k-1])
                            gamma = material.Γ
                            @inbounds ab, b5 = _diffusion_centralDifference_neighbors_(
                                gamma, mesh.dz[k], mesh.dz[k-1], (mesh.dx[i] * mesh.dy[j])
                            )
                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i,j,k-1]
                            AV[n] = ab
                            b[id] += b5
                        end

                        #Top Coefficents
                        if (k != mesh.n1) && (mesh.n1 != 1) && (phi.onoff[i,j,k+1])
                            gamma = material.Γ
                            @inbounds at, b6 = _diffusion_centralDifference_neighbors_(
                                gamma, mesh.dz[k], mesh.dz[k+1], (mesh.dx[i] * mesh.dy[j])
                            )
                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i,j,k+1]
                            AV[n] = at
                            b[id] += b6
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
                            n += 1
                            AI[n] = id
                            AJ[n] = id
                            AV[n] = ac - (aw + ae + as + an + ab + at)
                            b[id] += b0
                        else
                            n += 1
                            AI[n] = id
                            AJ[n] = id
                            AV[n] = -1.0 * (aw + ae + as + an + ab + at)
                        end
                    end
                end
            end
        end
    elseif !sparrays
        for i in 1:mesh.l1
            for j in 1:mesh.m1
                for k in 1:mesh.n1
                    if phi.onoff[i,j,k]
                        id = phi.gIndex[i,j,k]

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
                        if (i != 1) && (mesh.l1 != 1) && (phi.onoff[i-1,j,k])
                            gamma = material.Γ
                            @inbounds aw, b1 = _diffusion_centralDifference_neighbors_(
                                gamma, mesh.dx[i], mesh.dx[i-1], (mesh.dy[j] * mesh.dz[k])
                            )
                            A[id, phi.gIndex[i-1,j,k]] = aw
                            b[id] += b1
                        end

                        #East Coefficents
                        if (i != mesh.l1) && (mesh.l1 != 1) && (phi.onoff[i+1,j,k])
                            gamma = material.Γ
                            @inbounds ae, b2 = _diffusion_centralDifference_neighbors_(
                                gamma, mesh.dx[i], mesh.dx[i+1], (mesh.dy[j] * mesh.dz[k])
                            )
                            A[id, phi.gIndex[i+1,j,k]] = ae
                            b[id] += b2
                        end

                        #South Coefficents
                        if (j != 1) && (mesh.m1 != 1) && (phi.onoff[i,j-1,k])
                            gamma = material.Γ
                            @inbounds as, b3 = _diffusion_centralDifference_neighbors_(
                                gamma, mesh.dy[j], mesh.dy[j-1], (mesh.dx[i] * mesh.dz[k])
                            )
                            A[id, phi.gIndex[i,j-1,k]] = as
                            b[id] += b3
                        end

                        #North Coefficents
                        if (j != mesh.m1) && (mesh.m1 != 1) && (phi.onoff[i,j+1,k])
                            gamma = material.Γ
                            @inbounds an, b4 = _diffusion_centralDifference_neighbors_(
                                gamma, mesh.dy[j], mesh.dy[j+1], (mesh.dx[i] * mesh.dz[k])
                            )
                            A[id, phi.gIndex[i,j+1,k]] = an
                            b[id] += b4
                        end

                        #Bottom Coefficents
                        if (k != 1) && (mesh.n1 != 1) && (phi.onoff[i,j,k-1])
                            gamma = material.Γ
                            @inbounds ab, b5 = _diffusion_centralDifference_neighbors_(
                                gamma, mesh.dz[k], mesh.dz[k-1], (mesh.dx[i] * mesh.dy[j])
                            )
                            A[id, phi.gIndex[i,j,k-1]] = ab
                            b[id] += b5
                        end

                        #Top Coefficents
                        if (k != mesh.n1) && (mesh.n1 != 1) && (phi.onoff[i,j,k+1])
                            gamma = material.Γ
                            @inbounds at, b6 = _diffusion_centralDifference_neighbors_(
                                gamma, mesh.dz[k], mesh.dz[k+1], (mesh.dx[i] * mesh.dy[j])
                            )
                            A[id, phi.gIndex[i,j,k+1]] = at
                            b[id] += b6
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
                            A[id, id] = ac - (aw + ae + as + an + ab + at)
                            b[id] += b0
                        else
                            A[id, id] = -1.0 * (aw + ae + as + an + ab + at)
                        end
                    end
                end
            end
        end
    end

    if sparrays
        A = sparse(AI[1:n], AJ[1:n], AV[1:n])
    end

    return A, b
end
