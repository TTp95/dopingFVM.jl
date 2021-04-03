"""

"""
function _discretize_convection_downwind_ end

function _discretize_convection_downwind_(
    vel::CSVelocity1D,
    phi::CSPhi1D,
    bounds::Dict{String,BoundsStructured},
    material::CSMaterial1D,
    mesh::UnionCSMesh1D;
    velocityU::Array{<:AbstractFloat,1} = vel.fValues.uFace,
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
    mthreads::Bool = false,
    sparrays::Bool = true,
    interpolation::Signed = 1,
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
                awc = 0.0
                b1 = 0.0
                ae = 0.0
                aec = 0.0
                b2 = 0.0
                mflux = 0.0

                #West Coefficents
                if (i != 1) && (mesh.l1 != 1) && (phi.onoff[i-1])
                    rho = density_interpolation(
                        mesh.dx[i], mesh.dx[i-1], material.ρ[i], material.ρ[i-1];
                        interpolation = interpolation
                    )

                    mflux = -1.0 * rho * velocityU[i] * (1.0)
                    aw, awc, b1 = _convection_downwind_neighbors_(mflux)

                    n += 1
                    AI[n] = id
                    AJ[n] = phi.gIndex[i-1]
                    AV[n] = aw
                    b[id] += b1
                end

                #East Coefficents
                if (i != mesh.l1) && (mesh.l1 != 1)  && (phi.onoff[i+1])
                    rho = density_interpolation(
                        mesh.dx[i], mesh.dx[i+1], material.ρ[i], material.ρ[i+1];
                        interpolation = interpolation
                    )

                    mflux = rho * velocityU[i+1] * (1.0)
                    ae, aec, b2 = _convection_downwind_neighbors_(mflux)

                    n += 1
                    AI[n] = id
                    AJ[n] = phi.gIndex[i+1]
                    AV[n] = ae
                    b[id] += b2
                end

                #Center Coefficent
                if phi.bounds[i]
                    @inbounds ac, b0 = _convection_downwind_central_(
                        i,
                        velocityU,
                        phi,
                        bounds,
                        mesh;
                        T = T,
                    )
                    n += 1
                    AI[n] = id
                    AJ[n] = id
                    AV[n] = ac + aec + awc
                    b[id] += b0
                else
                    n += 1
                    AI[n] = id
                    AJ[n] = id
                    AV[n] = ac + aec + awc
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
                awc = 0.0
                b1 = 0.0
                ae = 0.0
                aec = 0.0
                b2 = 0.0
                mflux = 0.0

                #West Coefficents
                if (i != 1) && (mesh.l1 != 1) && (phi.onoff[i-1])
                    rho = density_interpolation(
                        mesh.dx[i], mesh.dx[i-1], material.ρ[i], material.ρ[i-1];
                        interpolation = interpolation
                    )

                    mflux = -1.0 * rho * velocityU[i] * (1.0)
                    aw, awc, b1 = _convection_downwind_neighbors_(mflux)

                    A[id, phi.gIndex[i-1]] = aw
                    b[id] += b1
                end

                #East Coefficents
                if (i != mesh.l1) && (mesh.l1 != 1)  && (phi.onoff[i+1])
                    rho = density_interpolation(
                        mesh.dx[i], mesh.dx[i+1], material.ρ[i], material.ρ[i+1];
                        interpolation = interpolation
                    )

                    mflux = rho * velocityU[i+1] * (1.0)
                    ae, aec, b2 = _convection_downwind_neighbors_(mflux)

                    A[id, phi.gIndex[i+1]] = ae
                    b[id] += b2
                end

                #Center Coefficent
                if phi.bounds[i]
                    @inbounds ac, b0 = _convection_downwind_central_(
                        i,
                        velocityU,
                        phi,
                        bounds,
                        mesh;
                        T = T,
                    )
                    A[id, id] = ac + aec + awc
                    b[id] += b0
                else
                    A[id, id] = ac + aec + awc
                end
            end
        end
    end

    if sparrays
        A = sparse(AI[1:n], AJ[1:n], AV[1:n])
    end

    return A, b
end

function _discretize_convection_downwind_(
    vel::CSVelocity1D,
    phi::CSPhi1D,
    bounds::Dict{String,BoundsStructured},
    material::UnionCSConstantMaterial,
    mesh::UnionCSMesh1D;
    velocityU::Array{<:AbstractFloat,1} = vel.fValues.uFace,
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
    mthreads::Bool = false,
    sparrays::Bool = true,
    interpolation::Signed = 1,
)
    n_equations = maximum_globalIndex(phi)

    if sparrays
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
                awc = 0.0
                b1 = 0.0
                ae = 0.0
                aec = 0.0
                b2 = 0.0
                mflux = 0.0

                #West Coefficents
                if (i != 1) && (mesh.l1 != 1) && (phi.onoff[i-1])
                    rho = material.ρ

                    mflux = -1.0 * rho * velocityU[i] * (1.0)
                    aw, awc, b1 = _convection_downwind_neighbors_(mflux)

                    n += 1
                    AI[n] = id
                    AJ[n] = phi.gIndex[i-1]
                    AV[n] = aw
                    b[id] += b1
                end

                #East Coefficents
                if (i != mesh.l1) && (mesh.l1 != 1)  && (phi.onoff[i+1])
                    rho = material.ρ

                    mflux = rho * velocityU[i+1] * (1.0)
                    ae, aec, b2 = _convection_downwind_neighbors_(mflux)

                    n += 1
                    AI[n] = id
                    AJ[n] = phi.gIndex[i+1]
                    AV[n] = ae
                    b[id] += b2
                end

                #Center Coefficent
                if phi.bounds[i]
                    @inbounds ac, b0 = _convection_downwind_central_(
                        i,
                        velocityU,
                        phi,
                        bounds,
                        mesh;
                        T = T,
                    )
                    n += 1
                    AI[n] = id
                    AJ[n] = id
                    AV[n] = ac + aec + awc
                    b[id] += b0
                else
                    n += 1
                    AI[n] = id
                    AJ[n] = id
                    AV[n] = ac + aec + awc
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
                awc = 0.0
                b1 = 0.0
                ae = 0.0
                aec = 0.0
                b2 = 0.0
                mflux = 0.0

                #West Coefficents
                if (i != 1) && (mesh.l1 != 1) && (phi.onoff[i-1])
                    rho = material.ρ

                    mflux = -1.0 * rho * velocityU[i] * (1.0)
                    aw, awc, b1 = _convection_downwind_neighbors_(mflux)

                    A[id, phi.gIndex[i-1]] = aw
                    b[id] += b1
                end

                #East Coefficents
                if (i != mesh.l1) && (mesh.l1 != 1)  && (phi.onoff[i+1])
                    rho = material.ρ

                    mflux = rho * velocityU[i+1] * (1.0)
                    ae, aec, b2 = _convection_downwind_neighbors_(mflux)

                    A[id, phi.gIndex[i+1]] = ae
                    b[id] += b2
                end

                #Center Coefficent
                if phi.bounds[i]
                    @inbounds ac, b0 = _convection_downwind_central_(
                        i,
                        velocityU,
                        phi,
                        bounds,
                        mesh;
                        T = T,
                    )
                    A[id, id] = ac + aec + awc
                    b[id] += b0
                else
                    A[id, id] = ac + aec + awc
                end
            end
        end
    end

    if sparrays
        A = sparse(AI[1:n], AJ[1:n], AV[1:n])
    end

    return A, b
end

function _discretize_convection_downwind_(
    vel::CSVelocity2D,
    phi::CSPhi2D,
    bounds::Dict{String,BoundsStructured},
    material::CSMaterial2D,
    mesh::UnionCSMesh2D;
    velocityU::Array{<:AbstractFloat,2} = vel.fValues.uFace,
    velocityV::Array{<:AbstractFloat,2} = vel.fValues.vFace,
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
    mthreads::Bool = false,
    sparrays::Bool = true,
    interpolation::Signed = 1,
)
    n_equations = maximum_globalIndex(phi)

    if sparrays
        n = 0
        AI = zeros(N, (5 * n_equations))
        AJ = zeros(N, (5 * n_equations))
        AV = zeros(T, (5 * n_equations))
    elseif !sparrays
        A = zeros(T, n_equations, n_equations)
    end

    b = zeros(T, n_equations)

    if mthreads
        #Base.Threads.@threads
    elseif sparrays
        for i in 1:mesh.l1
            for j in 1:mesh.m1
                if phi.onoff[i,j]
                    id = phi.gIndex[i,j]

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
                    mflux = 0.0

                    #West Coefficents
                    if (i != 1) && (mesh.l1 != 1) && (phi.onoff[i-1,j])
                        rho = density_interpolation(
                            mesh.dx[i], mesh.dx[i-1], material.ρ[i,j], material.ρ[i-1,j];
                            interpolation = interpolation
                        )

                        mflux = -1.0 * rho * velocityU[i,j] * (mesh.dy[j])
                        aw, awc, b1 = _convection_downwind_neighbors_(mflux)

                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i-1,j]
                        AV[n] = aw
                        b[id] += b1
                    end

                    #East Coefficents
                    if (i != mesh.l1) && (mesh.l1 != 1)  && (phi.onoff[i+1,j])
                        rho = density_interpolation(
                            mesh.dx[i], mesh.dx[i+1], material.ρ[i,j], material.ρ[i+1,j];
                            interpolation = interpolation
                        )

                        mflux = rho * velocityU[i+1,j] * (mesh.dy[j])
                        ae, aec, b2 = _convection_downwind_neighbors_(mflux)

                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i+1,j]
                        AV[n] = ae
                        b[id] += b2
                    end

                    #South Coefficents
                    if (j != 1) && (mesh.m1 != 1) && (phi.onoff[i,j-1])
                        rho = density_interpolation(
                            mesh.dy[j], mesh.dy[j-1], material.ρ[i,j], material.ρ[i,j-1];
                            interpolation = interpolation
                        )

                        mflux = -1.0 * rho * velocityV[i,j] * (mesh.dx[i])
                        as, asc, b3 = _convection_downwind_neighbors_(mflux)

                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i,j-1]
                        AV[n] = as
                        b[id] += b3
                    end

                    #North Coefficents
                    if (j != mesh.m1) && (mesh.m1 != 1) && (phi.onoff[i,j+1])
                        rho = density_interpolation(
                            mesh.dy[j], mesh.dy[j+1], material.ρ[i,j], material.ρ[i,j+1];
                            interpolation = interpolation
                        )

                        mflux = rho * velocityV[i,j+1] * (mesh.dx[i])
                        an, anc, b4 = _convection_downwind_neighbors_(mflux)

                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i,j+1]
                        AV[n] = an
                        b[id] += b4
                    end

                    #Center Coefficent
                    if phi.bounds[i,j]
                        @inbounds ac, b0 = _convection_downwind_central_(
                            i,
                            j,
                            velocityU,
                            velocityV,
                            phi,
                            bounds,
                            mesh;
                            T = T,
                        )
                        n += 1
                        AI[n] = id
                        AJ[n] = id
                        AV[n] = ac + (aec + awc + anc + asc)
                        b[id] += b0
                    else
                        n += 1
                        AI[n] = id
                        AJ[n] = id
                        AV[n] = ac + (aec + awc + anc + asc)
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
                    mflux = 0.0

                    #West Coefficents
                    if (i != 1) && (mesh.l1 != 1) && (phi.onoff[i-1,j])
                        rho = density_interpolation(
                            mesh.dx[i], mesh.dx[i-1], material.ρ[i,j], material.ρ[i-1,j];
                            interpolation = interpolation
                        )

                        mflux = -1.0 * rho * velocityU[i,j] * (mesh.dy[j])
                        aw, awc, b1 = _convection_downwind_neighbors_(mflux)

                        A[id, phi.gIndex[i-1,j]] = aw
                        b[id] += b1
                    end

                    #East Coefficents
                    if (i != mesh.l1) && (mesh.l1 != 1)  && (phi.onoff[i+1,j])
                        rho = density_interpolation(
                            mesh.dx[i], mesh.dx[i+1], material.ρ[i,j], material.ρ[i+1,j];
                            interpolation = interpolation
                        )

                        mflux = rho * velocityU[i+1,j] * (mesh.dy[j])
                        ae, aec, b2 = _convection_downwind_neighbors_(mflux)

                        A[id, phi.gIndex[i+1,j]] = ae
                        b[id] += b2
                    end

                    #South Coefficents
                    if (j != 1) && (mesh.m1 != 1) && (phi.onoff[i,j-1])
                        rho = density_interpolation(
                            mesh.dy[j], mesh.dy[j-1], material.ρ[i,j], material.ρ[i,j-1];
                            interpolation = interpolation
                        )

                        mflux = -1.0 * rho * velocityV[i,j] * (mesh.dx[i])
                        as, asc, b3 = _convection_downwind_neighbors_(mflux)

                        A[id, phi.gIndex[i,j-1]] = as
                        b[id] += b3
                    end

                    #North Coefficents
                    if (j != mesh.m1) && (mesh.m1 != 1) && (phi.onoff[i,j+1])
                        rho = density_interpolation(
                            mesh.dy[j], mesh.dy[j+1], material.ρ[i,j], material.ρ[i,j+1];
                            interpolation = interpolation
                        )

                        mflux = rho * velocityV[i,j+1] * (mesh.dx[i])
                        an, anc, b4 = _convection_downwind_neighbors_(mflux)

                        A[id, phi.gIndex[i,j+1]] = an
                        b[id] += b4
                    end

                    #Center Coefficent
                    if phi.bounds[i,j]
                        @inbounds ac, b0 = _convection_downwind_central_(
                            i,
                            j,
                            velocityU,
                            velocityV,
                            phi,
                            bounds,
                            mesh;
                            T = T,
                        )
                        A[id, id] = ac + aec + awc + anc + asc
                        b[id] += b0
                    else
                        A[id, id] = ac + aec + awc + anc + asc
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


function _discretize_convection_downwind_(
    vel::CSVelocity2D,
    phi::CSPhi2D,
    bounds::Dict{String,BoundsStructured},
    material::UnionCSConstantMaterial,
    mesh::UnionCSMesh2D;
    velocityU::Array{<:AbstractFloat,2} = vel.fValues.uFace,
    velocityV::Array{<:AbstractFloat,2} = vel.fValues.vFace,
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
    mthreads::Bool = false,
    sparrays::Bool = true,
    interpolation::Signed = 1,
)
    n_equations = maximum_globalIndex(phi)

    if sparrays
        n = 0
        AI = zeros(N, (5 * n_equations))
        AJ = zeros(N, (5 * n_equations))
        AV = zeros(T, (5 * n_equations))
    elseif !sparrays
        A = zeros(T, n_equations, n_equations)
    end

    b = zeros(T, n_equations)

    if mthreads
        #Base.Threads.@threads
    elseif sparrays
        for i in 1:mesh.l1
            for j in 1:mesh.m1
                if phi.onoff[i,j]
                    id = phi.gIndex[i,j]

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
                    mflux = 0.0

                    #West Coefficents
                    if (i != 1) && (mesh.l1 != 1) && (phi.onoff[i-1,j])
                        rho = material.ρ

                        mflux = -1.0 * rho * velocityU[i,j] * (mesh.dy[j])
                        aw, awc, b1 = _convection_downwind_neighbors_(mflux)

                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i-1,j]
                        AV[n] = aw
                        b[id] += b1
                    end

                    #East Coefficents
                    if (i != mesh.l1) && (mesh.l1 != 1)  && (phi.onoff[i+1,j])
                        rho = material.ρ

                        mflux = rho * velocityU[i+1,j] * (mesh.dy[j])
                        ae, aec, b2 = _convection_downwind_neighbors_(mflux)

                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i+1,j]
                        AV[n] = ae
                        b[id] += b2
                    end

                    #South Coefficents
                    if (j != 1) && (mesh.m1 != 1) && (phi.onoff[i,j-1])
                        rho = material.ρ

                        mflux = -1.0 * rho * velocityV[i,j] * (mesh.dx[i])
                        as, asc, b3 = _convection_downwind_neighbors_(mflux)

                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i,j-1]
                        AV[n] = as
                        b[id] += b3
                    end

                    #North Coefficents
                    if (j != mesh.m1) && (mesh.m1 != 1) && (phi.onoff[i,j+1])
                        rho = material.ρ

                        mflux = rho * velocityV[i,j+1] * (mesh.dx[i])
                        an, anc, b4 = _convection_downwind_neighbors_(mflux)

                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i,j+1]
                        AV[n] = an
                        b[id] += b4
                    end

                    #Center Coefficent
                    if phi.bounds[i,j]
                        @inbounds ac, b0 = _convection_downwind_central_(
                            i,
                            j,
                            velocityU,
                            velocityV,
                            phi,
                            bounds,
                            mesh;
                            T = T,
                        )
                        n += 1
                        AI[n] = id
                        AJ[n] = id
                        AV[n] = ac + aec + awc + anc + asc
                        b[id] += b0
                    else
                        n += 1
                        AI[n] = id
                        AJ[n] = id
                        AV[n] = ac + aec + awc + anc + asc
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
                    mflux = 0.0

                    #West Coefficents
                    if (i != 1) && (mesh.l1 != 1) && (phi.onoff[i-1,j])
                        rho = material.ρ

                        mflux = -1.0 * rho * velocityU[i,j] * (mesh.dy[j])
                        aw, awc, b1 = _convection_downwind_neighbors_(mflux)

                        A[id, phi.gIndex[i-1,j]] = aw
                        b[id] += b1
                    end

                    #East Coefficents
                    if (i != mesh.l1) && (mesh.l1 != 1)  && (phi.onoff[i+1,j])
                        rho = material.ρ

                        mflux = rho * velocityU[i+1,j] * (mesh.dy[j])
                        ae, aec, b2 = _convection_downwind_neighbors_(mflux)

                        A[id, phi.gIndex[i+1,j]] = ae
                        b[id] += b2
                    end

                    #South Coefficents
                    if (j != 1) && (mesh.m1 != 1) && (phi.onoff[i,j-1])
                        rho = material.ρ

                        mflux = -1.0 * rho * velocityV[i,j] * (mesh.dx[i])
                        as, asc, b3 = _convection_downwind_neighbors_(mflux)

                        A[id, phi.gIndex[i,j-1]] = as
                        b += b3
                    end

                    #North Coefficents
                    if (j != mesh.m1) && (mesh.m1 != 1) && (phi.onoff[i,j+1])
                        rho = material.ρ

                        mflux = rho * velocityV[i,j+1] * (mesh.dx[i])
                        an, anc, b4 = _convection_downwind_neighbors_(mflux)

                        A[id, phi.gIndex[i,j+1]] = an
                        b[id] += b4
                    end

                    #Center Coefficent
                    if phi.bounds[i,j]
                        @inbounds ac, b0 = _convection_downwind_central_(
                            i,
                            j,
                            velocityU,
                            velocityV,
                            phi,
                            bounds,
                            mesh;
                            T = T,
                        )
                        A[id, id] = ac + aec + awc + anc + asc
                        b[id] += b0
                    else
                        A[id, id] = ac + aec + awc + anc + asc
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

function _discretize_convection_downwind_(
    vel::CSVelocity3D,
    phi::CSPhi3D,
    bounds::Dict{String,BoundsStructured},
    material::CSMaterial3D,
    mesh::UnionCSMesh3D;
    velocityU::Array{<:AbstractFloat,3} = vel.fValues.uFace,
    velocityV::Array{<:AbstractFloat,3} = vel.fValues.vFace,
    velocityW::Array{<:AbstractFloat,3} = vel.fValues.wFace,
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
    mthreads::Bool = false,
    sparrays::Bool = true,
    interpolation::Signed = 1,
)
    n_equations = maximum_globalIndex(phi)

    if sparrays
        n = 0
        AI = zeros(N, (7 * n_equations))
        AJ = zeros(N, (7 * n_equations))
        AV = zeros(T, (7 * n_equations))
    elseif !sparrays
        A = zeros(T, n_equations, n_equations)
    end

    b = zeros(T, n_equations)

    if mthreads
        #Base.Threads.@threads
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
                        mflux = 0.0

                        #West Coefficents
                        if (i != 1) && (mesh.l1 != 1) && (phi.onoff[i-1,j,k])
                            rho = density_interpolation(
                                mesh.dx[i], mesh.dx[i-1], material.ρ[i,j,k], material.ρ[i-1,j,k];
                                interpolation = interpolation
                            )

                            mflux = -1.0 * rho * velocityU[i,j,k] * (mesh.dy[j] * mesh.dz[k])
                            aw, awc, b1 = _convection_downwind_neighbors_(mflux)

                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i-1,j,k]
                            AV[n] = aw
                            b[id] += b1
                        end

                        #East Coefficents
                        if (i != mesh.l1) && (mesh.l1 != 1)  && (phi.onoff[i+1,j,k])
                            rho = density_interpolation(
                                mesh.dx[i], mesh.dx[i+1], material.ρ[i,j,k], material.ρ[i+1,j,k];
                                interpolation = interpolation
                            )

                            mflux = rho * velocityU[i+1,j,k] * (mesh.dy[j] * mesh.dz[k])
                            ae, aec, b2 = _convection_downwind_neighbors_(mflux)

                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i+1,j,k]
                            AV[n] = ae
                            b[id] += b2
                        end

                        #South Coefficents
                        if (j != 1) && (mesh.m1 != 1) && (phi.onoff[i,j-1,k])
                            rho = density_interpolation(
                                mesh.dy[j], mesh.dy[j-1], material.ρ[i,j,k], material.ρ[i,j-1,k];
                                interpolation = interpolation
                            )

                            mflux = -1.0 * rho * velocityV[i,j,k] * (mesh.dx[i] * mesh.dz[k])
                            as, asc, b3 = _convection_downwind_neighbors_(mflux)

                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i,j-1,k]
                            AV[n] = as
                            b[id] += b3
                        end

                        #North Coefficents
                        if (j != mesh.m1) && (mesh.m1 != 1) && (phi.onoff[i,j+1,k])
                            rho = density_interpolation(
                                mesh.dy[j], mesh.dy[j+1], material.ρ[i,j,k], material.ρ[i,j+1,k];
                                interpolation = interpolation
                            )

                            mflux = rho * velocityV[i,j+1,k] * (mesh.dx[i] * mesh.dz[k])
                            an, anc, b4 = _convection_downwind_neighbors_(mflux)

                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i,j+1,k]
                            AV[n] = an
                            b[id] += b4
                        end

                        #Bottom Coefficents
                        if (k != 1) && (mesh.n1 != 1) && (phi.onoff[i,j,k-1])
                            rho = density_interpolation(
                                mesh.dz[k], mesh.dz[k-1], material.ρ[i,j,k], material.ρ[i,j,k-1];
                                interpolation = interpolation
                            )

                            mflux = -1.0 * rho * velocityW[i,j,k] * (mesh.dx[i] * mesh.dy[j])
                            ab, abc, b5 = _convection_downwind_neighbors_(mflux)

                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i,j,k-1]
                            AV[n] = ab
                            b[id] += b5
                        end

                        #Top Coefficents
                        if (k != mesh.n1) && (mesh.n1 != 1) && (phi.onoff[i,j,k+1])
                            rho = density_interpolation(
                                mesh.dz[k], mesh.dz[k+1], material.ρ[i,j,k], material.ρ[i,j,k+1];
                                interpolation = interpolation
                            )

                            mflux = rho * velocityW[i,j,k+1] * (mesh.dx[i] * mesh.dy[j])
                            at, atc, b6 = _convection_downwind_neighbors_(mflux)

                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i,j,k+1]
                            AV[n] = at
                            b[id] += b6
                        end

                        #Center Coefficent
                        if phi.bounds[i,j,k]
                            @inbounds ac, b0 = _convection_downwind_central_(
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
                            n += 1
                            AI[n] = id
                            AJ[n] = id
                            AV[n] = ac + aec + awc + anc + asc + atc + abc
                            b[id] += b0
                        else
                            n += 1
                            AI[n] = id
                            AJ[n] = id
                            AV[n] = ac + aec + awc + anc + asc + atc + abc
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
                        mflux = 0.0

                        #West Coefficents
                        if (i != 1) && (mesh.l1 != 1) && (phi.onoff[i-1,j,k])
                            rho = density_interpolation(
                                mesh.dx[i], mesh.dx[i-1], material.ρ[i,j,k], material.ρ[i-1,j,k];
                                interpolation = interpolation
                            )

                            mflux = -1.0 * rho * velocityU[i,j,k] * (mesh.dy[j] * mesh.dz[k])
                            aw, awc, b1 = _convection_downwind_neighbors_(mflux)

                            A[id, phi.gIndex[i-1,j,k]] = aw
                            b[id] += b1
                        end

                        #East Coefficents
                        if (i != mesh.l1) && (mesh.l1 != 1)  && (phi.onoff[i+1,j,k])
                            rho = density_interpolation(
                                mesh.dx[i], mesh.dx[i+1], material.ρ[i,j,k], material.ρ[i+1,j,k];
                                interpolation = interpolation
                            )

                            mflux = rho * velocityU[i+1,j,k] * (mesh.dy[j] * mesh.dz[k])
                            ae, aec, b2 = _convection_downwind_neighbors_(mflux)

                            A[id, phi.gIndex[i+1,j,k]] = ae
                            b[id] += b2
                        end

                        #South Coefficents
                        if (j != 1) && (mesh.m1 != 1) && (phi.onoff[i,j-1,k])
                            rho = density_interpolation(
                                mesh.dy[j], mesh.dy[j-1], material.ρ[i,j,k], material.ρ[i,j-1,k];
                                interpolation = interpolation
                            )

                            mflux = -1.0 * rho * velocityV[i,j,k] * (mesh.dx[i] * mesh.dz[k])
                            as, asc, b3 = _convection_downwind_neighbors_(mflux)

                            A[id, phi.gIndex[i,j-1,k]] = as
                            b[id] += b3
                        end

                        #North Coefficents
                        if (j != mesh.m1) && (mesh.m1 != 1) && (phi.onoff[i,j+1,k])
                            rho = density_interpolation(
                                mesh.dy[j], mesh.dy[j+1], material.ρ[i,j,k], material.ρ[i,j+1,k];
                                interpolation = interpolation
                            )

                            mflux = rho * velocityV[i,j+1,k] * (mesh.dx[i] * mesh.dz[k])
                            an, anc, b4 = _convection_downwind_neighbors_(mflux)

                            A[id, phi.gIndex[i,j+1,k]] = an
                            b[id] += b4
                        end

                        #Bottom Coefficents
                        if (k != 1) && (mesh.n1 != 1) && (phi.onoff[i,j,k-1])
                            rho = density_interpolation(
                                mesh.dz[k], mesh.dz[k-1], material.ρ[i,j,k], material.ρ[i,j,k-1];
                                interpolation = interpolation
                            )

                            mflux = -1.0 * rho * velocityW[i,j,k] * (mesh.dx[i] * mesh.dy[j])
                            ab, abc, b5 = _convection_downwind_neighbors_(mflux)

                            A[id, phi.gIndex[i,j,k-1]] = ab
                            b[id] += b5
                        end

                        #Top Coefficents
                        if (k != mesh.n1) && (mesh.n1 != 1) && (phi.onoff[i,j,k+1])
                            rho = density_interpolation(
                                mesh.dz[k], mesh.dz[k+1], material.ρ[i,j,k], material.ρ[i,j,k+1];
                                interpolation = interpolation
                            )

                            mflux = rho * velocityW[i,j,k+1] * (mesh.dx[i] * mesh.dy[j])
                            at, atc, b6 = _convection_downwind_neighbors_(mflux)

                            A[id, phi.gIndex[i,j,k+1]] = at
                            b[id] += b6
                        end

                        #Center Coefficent
                        if phi.bounds[i,j,k]
                            @inbounds ac, b0 = _convection_downwind_central_(
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
                            A[id, id] = ac + aec + awc + anc + asc + atc + abc
                            b[id] += b0
                        else
                            A[id, id] = ac + aec + awc + anc + asc + atc + abc
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

function _discretize_convection_downwind_(
    vel::CSVelocity3D,
    phi::CSPhi3D,
    bounds::Dict{String,BoundsStructured},
    material::UnionCSConstantMaterial,
    mesh::UnionCSMesh3D;
    velocityU::Array{<:AbstractFloat,3} = vel.fValues.uFace,
    velocityV::Array{<:AbstractFloat,3} = vel.fValues.vFace,
    velocityW::Array{<:AbstractFloat,3} = vel.fValues.wFace,
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
    mthreads::Bool = false,
    sparrays::Bool = true,
    interpolation::Signed = 1,
)
    n_equations = maximum_globalIndex(phi)

    if sparrays
        n = 0
        AI = zeros(N, (7 * n_equations))
        AJ = zeros(N, (7 * n_equations))
        AV = zeros(T, (7 * n_equations))
    elseif !sparrays
        A = zeros(T, n_equations, n_equations)
    end

    b = zeros(T, n_equations)

    if mthreads
        #Base.Threads.@threads
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
                        mflux = 0.0

                        #West Coefficents
                        if (i != 1) && (mesh.l1 != 1) && (phi.onoff[i-1,j,k])
                            rho = material.ρ

                            mflux = -1.0 * rho * velocityU[i,j,k] * (mesh.dy[j] * mesh.dz[k])
                            aw, awc, b1 = _convection_downwind_neighbors_(mflux)

                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i-1,j,k]
                            AV[n] = aw
                            b[id] += b1
                        end

                        #East Coefficents
                        if (i != mesh.l1) && (mesh.l1 != 1)  && (phi.onoff[i+1,j,k])
                            rho = material.ρ

                            mflux = rho * velocityU[i+1,j,k] * (mesh.dy[j] * mesh.dz[k])
                            ae, aec, b2 = _convection_downwind_neighbors_(mflux)

                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i+1,j,k]
                            AV[n] = ae
                            b[id] += b2
                        end

                        #South Coefficents
                        if (j != 1) && (mesh.m1 != 1) && (phi.onoff[i,j-1,k])
                            rho = material.ρ

                            mflux = -1.0 * rho * velocityV[i,j,k] * (mesh.dx[i] * mesh.dz[k])
                            as, asc, b3 = _convection_downwind_neighbors_(mflux)

                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i,j-1,k]
                            AV[n] = as
                            b[id] += b3
                        end

                        #North Coefficents
                        if (j != mesh.m1) && (mesh.m1 != 1) && (phi.onoff[i,j+1,k])
                            rho = material.ρ

                            mflux = rho * velocityV[i,j+1,k] * (mesh.dx[i] * mesh.dz[k])
                            an, anc, b4 = _convection_downwind_neighbors_(mflux)

                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i,j+1,k]
                            AV[n] = an
                            b[id] += b4
                        end

                        #Bottom Coefficents
                        if (k != 1) && (mesh.n1 != 1) && (phi.onoff[i,j,k-1])
                            rho = material.ρ

                            mflux = -1.0 * rho * velocityW[i,j,k] * (mesh.dx[i] * mesh.dy[j])
                            ab, abc, b5 = _convection_downwind_neighbors_(mflux)

                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i,j,k-1]
                            AV[n] = ab
                            b[id] += b5
                        end

                        #Top Coefficents
                        if (k != mesh.n1) && (mesh.n1 != 1) && (phi.onoff[i,j,k+1])
                            rho = material.ρ

                            mflux = rho * velocityW[i,j,k+1] * (mesh.dx[i] * mesh.dy[j])
                            at, atc, b6 = _convection_downwind_neighbors_(mflux)

                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i,j,k+1]
                            AV[n] = at
                            b[id] += b6
                        end

                        #Center Coefficent
                        if phi.bounds[i,j,k]
                            @inbounds ac, b0 = _convection_downwind_central_(
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
                            n += 1
                            AI[n] = id
                            AJ[n] = id
                            AV[n] = ac + aec + awc + anc + asc + atc + abc
                            b[id] += b0
                        else
                            n += 1
                            AI[n] = id
                            AJ[n] = id
                            AV[n] = ac + aec + awc + anc + asc + atc + abc
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
                        mflux = 0.0

                        #West Coefficents
                        if (i != 1) && (mesh.l1 != 1) && (phi.onoff[i-1,j,k])
                            rho = material.ρ

                            mflux = -1.0 * rho * velocityU[i,j,k] * (mesh.dy[j] * mesh.dz[k])
                            aw, awc, b1 = _convection_downwind_neighbors_(mflux)

                            A[id, phi.gIndex[i-1,j,k]] = aw
                            b[id] += b1
                        end

                        #East Coefficents
                        if (i != mesh.l1) && (mesh.l1 != 1)  && (phi.onoff[i+1,j,k])
                            rho = material.ρ

                            mflux = rho * velocityU[i+1,j,k] * (mesh.dy[j] * mesh.dz[k])
                            ae, aec, b2 = _convection_downwind_neighbors_(mflux)

                            A[id, phi.gIndex[i+1,j,k]] = ae
                            b[id] += b2
                        end

                        #South Coefficents
                        if (j != 1) && (mesh.m1 != 1) && (phi.onoff[i,j-1,k])
                            rho = material.ρ

                            mflux = -1.0 * rho * velocityV[i,j,k] * (mesh.dx[i] * mesh.dz[k])
                            as, asc, b3 = _convection_downwind_neighbors_(mflux)

                            A[id, phi.gIndex[i,j-1,k]] = as
                            b[id] += b3
                        end

                        #North Coefficents
                        if (j != mesh.m1) && (mesh.m1 != 1) && (phi.onoff[i,j+1,k])
                            rho = material.ρ

                            mflux = rho * velocityV[i,j+1,k] * (mesh.dx[i] * mesh.dz[k])
                            an, anc, b4 = _convection_downwind_neighbors_(mflux)

                            A[id, phi.gIndex[i,j+1,k]] = an
                            b[id] += b4
                        end

                        #Bottom Coefficents
                        if (k != 1) && (mesh.n1 != 1) && (phi.onoff[i,j,k-1])
                            rho = material.ρ

                            mflux = -1.0 * rho * velocityW[i,j,k] * (mesh.dx[i] * mesh.dy[j])
                            ab, abc, b5 = _convection_downwind_neighbors_(mflux)

                            A[id, phi.gIndex[i,j,k-1]] = ab
                            b[id] += b5
                        end

                        #Top Coefficents
                        if (k != mesh.n1) && (mesh.n1 != 1) && (phi.onoff[i,j,k+1])
                            rho = material.ρ

                            mflux = rho * velocityW[i,j,k+1] * (mesh.dx[i] * mesh.dy[j])
                            at, atc, b6 = _convection_downwind_neighbors_(mflux)

                            A[id, phi.gIndex[i,j,k+1]] = at
                            b[id] += b6
                        end

                        #Center Coefficent
                        if phi.bounds[i,j,k]
                            @inbounds ac, b0 = _convection_downwind_central_(
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
                            A[id, id] = ac + aec + awc + anc + asc + atc + abc
                            b[id] += b0
                        else
                            A[id, phi.gIndex[i,j,k]] = ac + aec + awc + anc + asc + atc + abc
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
