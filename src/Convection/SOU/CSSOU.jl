"""

"""
function _discretize_convection_secondorderupwind_ end

function _discretize_convection_secondorderupwind_(
    vel::CSVelocity1D,
    phi::CSPhi1D,
    bounds::Dict{String,BoundsStructured},
    material::CSMaterial1D,
    mesh::UnionCSMesh1D,
    inout::Bool = false;
    velocityU::Array{<:AbstractFloat,1} = vel.fValues.uFace,
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
    interpolation::Signed = 1,
)
    n_equations = maximum_globalIndex(phi)

    n = 0
    AI = zeros(N, (5 * n_equations))
    AJ = zeros(N, (5 * n_equations))
    AV = zeros(T, (5 * n_equations))

    b = zeros(T, n_equations)

    for i in 1:mesh.l1
        if phi.onoff[i]
            id = phi.gIndex[i]

            #Auxiliar variables
            ac = 0.0
            b0 = 0.0
            aw = 0.0
            awe = 0.0
            aww = 0.0
            awc = 0.0
            b1 = 0.0
            ae = 0.0
            aew = 0.0
            aee = 0.0
            aec = 0.0
            b2 = 0.0
            mflux = 0.0
            num = 0.0
            den = 0.0

            #West Coefficents
            if (i != 1) && (mesh.l1 != 1) && (phi.onoff[i-1])
                rho = density_interpolation(
                    mesh.dx[i], mesh.dx[i-1], material.ρ[i], material.ρ[i-1];
                    interpolation = interpolation
                )

                mflux = -1.0 * rho * velocityU[i] * (1.0)

                if (mflux >= 0) && (i != mesh.l1) && (phi.onoff[i+1])
                    num = ((mesh.x[i] - 0.5 * mesh.dx[i]) - mesh.x[i])
                    den = (mesh.x[i] - mesh.x[i+1])
                    alpha = abs(num / den)
                    awc, awe, b1 = _convection_secondorderupwind_neighbors_(mflux, alpha)
                    n += 1
                    AI[n] = id
                    AJ[n] = phi.gIndex[i+1]
                    AV[n] = awe
                    b[id] += b1

                elseif (mflux >= 0) && (i == mesh.l1)
                    num = ((mesh.x[i] - 0.5 * mesh.dx[i]) - mesh.x[i])
                    den = (mesh.x[i] - (mesh.x[i] + 0.5 * mesh.dx[i]))
                    alpha = abs(num / den)

                    if phi.bounds[i]
                        boundvalue = find_bondValue(i, id, phi, mesh, bounds, 'e')
                    else
                        boundvalue = phi.eval[i]
                    end

                    awc = (1.0 + alpha) * max(mflux, 0.0)
                    awe = 0.0
                    b1 = (alpha * boundvalue) * max(mflux, 0.0)

                    b[id] += b1

                elseif (mflux >= 0) && (!phi.onoff[i+1])
                    num = ((mesh.x[i] - 0.5 * mesh.dx[i]) - mesh.x[i])
                    den = (mesh.x[i] - (mesh.x[i] + 0.5 * mesh.dx[i]))
                    alpha = abs(num / den)

                    if phi.bounds[i]
                        boundvalue = find_bondValue(i, id, phi, mesh, bounds, 'e')
                    else
                        boundvalue = phi.eval[i+1]
                    end

                    awc = (1.0 + alpha) * max(mflux, 0.0)
                    awe = 0.0
                    b1 = (alpha * boundvalue) * max(mflux, 0.0)

                    b[id] += b1

                end


                if (mflux < 0) && (i != 2) && (phi.onoff[i-2])
                    num = ((mesh.x[i] - 0.5 * mesh.dx[i]) - mesh.x[i-1])
                    den = (mesh.x[i-1] - mesh.x[i-2])
                    alpha = abs(num / den)
                    aux_mflux = -1.0 * mflux
                    aw, aww, b1 = _convection_secondorderupwind_neighbors_(aux_mflux, alpha)
                    n += 1
                    AI[n] = id
                    AJ[n] = phi.gIndex[i-1]
                    AV[n] = -1.0 * aw

                    n += 1
                    AI[n] = id
                    AJ[n] = phi.gIndex[i-2]
                    AV[n] = -1.0 * aww

                    b[id] += b1

                elseif (mflux < 0) && (i == 2)
                    num = ((mesh.x[i] - 0.5 * mesh.dx[i]) - mesh.x[i-1])
                    den = (mesh.x[i-1] - (mesh.x[i-1] - 0.5 * mesh.dx[i-1]))
                    alpha = abs(num / den)
                    aux_mflux = -1.0 * mflux

                    if phi.bounds[i-1]
                        boundvalue = find_bondValue((i-1), (phi.gIndex[i-1]), phi, mesh, bounds, 'w')
                    else
                        boundvalue = phi.eval[i-1]
                    end

                    aw = (1.0 + alpha) * max(aux_mflux, 0.0)
                    aww = 0.0
                    b1 = (alpha * boundvalue) * max(aux_mflux, 0.0)

                    n += 1
                    AI[n] = id
                    AJ[n] = phi.gIndex[i-1]
                    AV[n] = -1.0 * aw

                    b[id] += -1.0 * b1

                elseif (mflux < 0) && (!phi.onoff[i-2])
                    num = ((mesh.x[i] - 0.5 * mesh.dx[i]) - mesh.x[i-1])
                    den = (mesh.x[i-1] - (mesh.x[i-1] - 0.5 * mesh.dx[i-1]))
                    alpha = abs(num / den)
                    aux_mflux = -1.0 * mflux

                    if phi.bounds[i-1]
                        boundvalue = find_bondValue((i-1), (phi.gIndex[i-1]), phi, mesh, bounds, 'w')
                    else
                        boundvalue = phi.eval[i-2]
                    end

                    aw = (1.0 + alpha) * max(aux_mflux, 0.0)
                    aww = 0.0
                    b1 = (alpha * boundvalue) * max(aux_mflux, 0.0)

                    n += 1
                    AI[n] = id
                    AJ[n] = phi.gIndex[i-1]
                    AV[n] = -1.0 * aw

                    b[id] += -1.0 * b1

                end

            end

            #East Coefficents
            if (i != mesh.l1) && (mesh.l1 != 1) && (phi.onoff[i+1])
                rho = density_interpolation(
                    mesh.dx[i], mesh.dx[i+1], material.ρ[i], material.ρ[i+1];
                    interpolation = interpolation
                )

                mflux = rho * velocityU[i+1] * (1.0)

                if (mflux >= 0) && (i != 1) && (phi.onoff[i-1])
                    num = ((mesh.x[i] + 0.5 * mesh.dx[i]) - mesh.x[i])
                    den = (mesh.x[i] - mesh.x[i-1])
                    alpha = abs(num / den)
                    aec, aew, b2 = _convection_secondorderupwind_neighbors_(mflux, alpha)
                    n += 1
                    AI[n] = id
                    AJ[n] = phi.gIndex[i-1]
                    AV[n] = aew
                    b[id] += b2

                elseif (mflux >= 0) && (i == 1)
                    num = ((mesh.x[i] + 0.5 * mesh.dx[i]) - mesh.x[i])
                    den = (mesh.x[i] - (mesh.x[i] - 0.5 * mesh.dx[i]))
                    alpha = abs(num / den)

                    if phi.bounds[i]
                        boundvalue = find_bondValue(i, id, phi, mesh, bounds, 'w')
                    else
                        boundvalue = phi.eval[i]
                    end

                    aec = (1.0 + alpha) * max(mflux, 0.0)
                    aew = 0.0
                    b2 = (alpha * boundvalue) * max(mflux, 0.0)

                    b[id] += b2

                elseif (mflux >= 0) && (!phi.onoff[i-1])
                    num = ((mesh.x[i] + 0.5 * mesh.dx[i]) - mesh.x[i])
                    den = (mesh.x[i] - (mesh.x[i] - 0.5 * mesh.dx[i]))
                    alpha = abs(num / den)

                    if phi.bounds[i]
                        boundvalue = find_bondValue(i, id, phi, mesh, bounds, 'w')
                    else
                        boundvalue = phi.eval[i-1]
                    end

                    aec = (1.0 + alpha) * max(mflux, 0.0)
                    aew = 0.0
                    b2 = (alpha * boundvalue) * max(mflux, 0.0)

                    b[id] += b2

                end


                if (mflux < 0) && (i != (mesh.l1-1)) && (phi.onoff[i+2])
                    num = ((mesh.x[i] + 0.5 * mesh.dx[i]) - mesh.x[i+1])
                    den = (mesh.x[i+1] - mesh.x[i+2])
                    alpha = abs(num / den)
                    aux_mflux = -1.0 * mflux
                    ae, aee, b2 = _convection_secondorderupwind_neighbors_(aux_mflux, alpha)
                    n += 1
                    AI[n] = id
                    AJ[n] = phi.gIndex[i+1]
                    AV[n] = -1.0 * ae

                    n += 1
                    AI[n] = id
                    AJ[n] = phi.gIndex[i+2]
                    AV[n] = -1.0 * aee

                    b[id] += b2

                elseif (mflux < 0) && (i == (mesh.l1-1))
                    num = ((mesh.x[i] + 0.5 * mesh.dx[i]) - mesh.x[i+1])
                    den = (mesh.x[i+1] - (mesh.x[i+1] + 0.5 * mesh.dx[i+1]))
                    alpha = abs(num / den)
                    aux_mflux = -1.0 * mflux

                    if phi.bounds[i+1]
                        boundvalue = find_bondValue((i+1), (phi.gIndex[i+1]), phi, mesh, bounds, 'e')
                    else
                        boundvalue = phi.eval[i+1]
                    end

                    ae = (1.0 + alpha) * max(aux_mflux, 0.0)
                    aee = 0.0
                    b2 = (alpha * boundvalue) * max(aux_mflux, 0.0)

                    n += 1
                    AI[n] = id
                    AJ[n] = phi.gIndex[i+1]
                    AV[n] = -1.0 * ae

                    b[id] += -1.0 * b2

                elseif (mflux < 0) && (!phi.onoff[i+2])
                    num = ((mesh.x[i] + 0.5 * mesh.dx[i]) - mesh.x[i+1])
                    den = (mesh.x[i+1] - (mesh.x[i+1] + 0.5 * mesh.dx[i+1]))
                    alpha = abs(num / den)
                    aux_mflux = -1.0 * mflux

                    if phi.bounds[i+1]
                        boundvalue = find_bondValue((i+1), (phi.gIndex[i+1]), phi, mesh, bounds, 'e')
                    else
                        boundvalue = phi.eval[i+2]
                    end

                    ae = (1.0 + alpha) * max(aux_mflux, 0.0)
                    aee = 0.0
                    b2 = (alpha * boundvalue) * max(aux_mflux, 0.0)

                    n += 1
                    AI[n] = id
                    AJ[n] = phi.gIndex[i+1]
                    AV[n] = -1.0 * ae

                    b[id] += -1.0 * b2

                end

            end

            #Center Coefficent
            if phi.bounds[i]
                @inbounds ac, b0 = _convection_secondorderupwind_central_(
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
                if inout
                    AV[n] = ac + aec + awc
                else
                    AV[n] = ac + (ae + aee + aw + aww) - (aew + awe) #ac + aec + awc
                end
                b[id] += b0
            else
                n += 1
                AI[n] = id
                AJ[n] = id
                if inout
                    AV[n] = ac + aec + awc
                else
                    AV[n] = ac + (ae + aee + aw + aww) - (aew + awe) #ac + aec + awc
                end
            end
        end
    end

    A = sparse(AI[1:n], AJ[1:n], AV[1:n])

    return A, b
end

function _discretize_convection_secondorderupwind_(
    vel::CSVelocity2D,
    phi::CSPhi2D,
    bounds::Dict{String,BoundsStructured},
    material::CSMaterial2D,
    mesh::UnionCSMesh2D,
    inout::Bool = false;
    velocityU::Array{<:AbstractFloat,2} = vel.fValues.uFace,
    velocityV::Array{<:AbstractFloat,2} = vel.fValues.vFace,
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
    interpolation::Signed = 1,
)
    n_equations = maximum_globalIndex(phi)

    n = 0
    AI = zeros(N, (9 * n_equations))
    AJ = zeros(N, (9 * n_equations))
    AV = zeros(T, (9 * n_equations))

    b = zeros(T, n_equations)

    for i in 1:mesh.l1
        for j in 1:mesh.m1
            if phi.onoff[i,j]
                id = phi.gIndex[i,j]

                #Auxiliar variables
                ac = 0.0
                b0 = 0.0
                aw = 0.0
                awe = 0.0
                aww = 0.0
                awc = 0.0
                b1 = 0.0
                ae = 0.0
                aew = 0.0
                aee = 0.0
                aec = 0.0
                b2 = 0.0
                as = 0.0
                asn = 0.0
                ass = 0.0
                asc = 0.0
                b3 = 0.0
                an = 0.0
                ans = 0.0
                ann = 0.0
                anc = 0.0
                b4 = 0.0
                mflux = 0.0
                mfluxw = 0.0
                mfluxe = 0.0
                mfluxs = 0.0
                mfluxn = 0.0
                acfw = 0.0
                acfe = 0.0
                acfs = 0.0
                acfn = 0.0
                num = 0.0
                den = 0.0

                #West Coefficents
                if (i != 1) && (mesh.l1 != 1) && (phi.onoff[i-1,j])
                    rho = density_interpolation(
                        mesh.dx[i], mesh.dx[i-1], material.ρ[i,j], material.ρ[i-1,j];
                        interpolation = interpolation
                    )

                    mflux = -1.0 * rho * velocityU[i,j] * (mesh.dy[j])
                    mfluxw = -1.0 * rho * velocityU[i,j] * (mesh.dy[j])

                    if (mflux >= 0) && (i != mesh.l1) && (phi.onoff[i+1,j])
                        num = ((mesh.x[i] - 0.5 * mesh.dx[i]) - mesh.x[i])
                        den = (mesh.x[i] - mesh.x[i+1])
                        alpha = abs(num / den)
                        awc, awe, b1 = _convection_secondorderupwind_neighbors_(mflux, alpha)
                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i+1,j]
                        AV[n] = awe
                        b[id] += b1

                    elseif (mflux >= 0) && (i == mesh.l1)
                        num = ((mesh.x[i] - 0.5 * mesh.dx[i]) - mesh.x[i])
                        den = (mesh.x[i] - (mesh.x[i] + 0.5 * mesh.dx[i]))
                        alpha = abs(num / den)

                        if phi.bounds[i,j]
                            boundvalue = find_bondValue(i, j, id, phi, mesh, bounds, 'e')
                        else
                            boundvalue = phi.eval[i,j]
                        end

                        awc = (1.0 + alpha) * max(mflux, 0.0)
                        acfw = (1.0 + alpha) * max(mflux, 0.0)
                        awe = 0.0
                        b1 = (alpha * boundvalue) * max(mflux, 0.0)

                        b[id] += b1

                    elseif (mflux >= 0) && (!phi.onoff[i+1,j])
                        num = ((mesh.x[i] - 0.5 * mesh.dx[i]) - mesh.x[i])
                        den = (mesh.x[i] - (mesh.x[i] + 0.5 * mesh.dx[i]))
                        alpha = abs(num / den)

                        if phi.bounds[i,j]
                            boundvalue = find_bondValue(i, j, id, phi, mesh, bounds, 'e')
                        else
                            boundvalue = phi.eval[i+1,j]
                        end

                        awc = (1.0 + alpha) * max(mflux, 0.0)
                        acfw = (1.0 + alpha) * max(mflux, 0.0)
                        awe = 0.0
                        b1 = (alpha * boundvalue) * max(mflux, 0.0)

                        b[id] += b1

                    end


                    if (mflux < 0) && (i != 2) && (phi.onoff[i-2,j])
                        num = ((mesh.x[i] - 0.5 * mesh.dx[i]) - mesh.x[i-1])
                        den = (mesh.x[i-1] - mesh.x[i-2])
                        alpha = abs(num / den)
                        aux_mflux = -1.0 * mflux
                        aw, aww, b1 = _convection_secondorderupwind_neighbors_(aux_mflux, alpha)
                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i-1,j]
                        AV[n] = -1.0 * aw

                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i-2,j]
                        AV[n] = -1.0 * aww

                        b[id] += b1

                    elseif (mflux < 0) && (i == 2)
                        num = ((mesh.x[i] - 0.5 * mesh.dx[i]) - mesh.x[i-1])
                        den = (mesh.x[i-1] - (mesh.x[i-1] - 0.5 * mesh.dx[i-1]))
                        alpha = abs(num / den)
                        aux_mflux = -1.0 * mflux

                        if phi.bounds[i-1,j]
                            boundvalue = find_bondValue((i-1), j, (phi.gIndex[i-1,j]), phi, mesh, bounds, 'w')
                        else
                            boundvalue = phi.eval[i-1,j]
                        end

                        aw = (1.0 + alpha) * max(aux_mflux, 0.0)
                        aww = 0.0
                        b1 = (alpha * boundvalue) * max(aux_mflux, 0.0)

                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i-1,j]
                        AV[n] = -1.0 * aw

                        b[id] += -1.0 * b1

                    elseif (mflux < 0) && (!phi.onoff[i-2,j])
                        num = ((mesh.x[i] - 0.5 * mesh.dx[i]) - mesh.x[i-1])
                        den = (mesh.x[i-1] - (mesh.x[i-1] - 0.5 * mesh.dx[i-1]))
                        alpha = abs(num / den)
                        aux_mflux = -1.0 * mflux

                        if phi.bounds[i-1,j]
                            boundvalue = find_bondValue((i-1), j, (phi.gIndex[i-1,j]), phi, mesh, bounds, 'w')
                        else
                            boundvalue = phi.eval[i-2,j]
                        end

                        aw = (1.0 + alpha) * max(aux_mflux, 0.0)
                        aww = 0.0
                        b1 = (alpha * boundvalue) * max(aux_mflux, 0.0)

                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i-1,j]
                        AV[n] = -1.0 * aw

                        b[id] += -1.0 * b1

                    end

                end

                #East Coefficents
                if (i != mesh.l1) && (mesh.l1 != 1) && (phi.onoff[i+1,j])
                    rho = density_interpolation(
                        mesh.dx[i], mesh.dx[i+1], material.ρ[i,j], material.ρ[i+1,j];
                        interpolation = interpolation
                    )

                    mflux = rho * velocityU[i+1,j] * (mesh.dy[j])
                    mfluxe = rho * velocityU[i+1,j] * (mesh.dy[j])

                    if (mflux >= 0) && (i != 1) && (phi.onoff[i-1,j])
                        num = ((mesh.x[i] + 0.5 * mesh.dx[i]) - mesh.x[i])
                        den = (mesh.x[i] - mesh.x[i-1])
                        alpha = abs(num / den)
                        aec, aew, b2 = _convection_secondorderupwind_neighbors_(mflux, alpha)
                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i-1,j]
                        AV[n] = aew
                        b[id] += b2

                    elseif (mflux >= 0) && (i == 1)
                        num = ((mesh.x[i] + 0.5 * mesh.dx[i]) - mesh.x[i])
                        den = (mesh.x[i] - (mesh.x[i] - 0.5 * mesh.dx[i]))
                        alpha = abs(num / den)

                        if phi.bounds[i,j]
                            boundvalue = find_bondValue(i, j, id, phi, mesh, bounds, 'w')
                        else
                            boundvalue = phi.eval[i,j]
                        end

                        aec = (1.0 + alpha) * max(mflux, 0.0)
                        acfe = (1.0 + alpha) * max(mflux, 0.0)
                        aew = 0.0
                        b2 = (alpha * boundvalue) * max(mflux, 0.0)

                        b[id] += b2

                    elseif (mflux >= 0) && (!phi.onoff[i-1,j])
                        num = ((mesh.x[i] + 0.5 * mesh.dx[i]) - mesh.x[i])
                        den = (mesh.x[i] - (mesh.x[i] - 0.5 * mesh.dx[i]))
                        alpha = abs(num / den)

                        if phi.bounds[i,j]
                            boundvalue = find_bondValue(i, j, id, phi, mesh, bounds, 'w')
                        else
                            boundvalue = phi.eval[i-1,j]
                        end

                        aec = (1.0 + alpha) * max(mflux, 0.0)
                        acfe = (1.0 + alpha) * max(mflux, 0.0)
                        aew = 0.0
                        b2 = (alpha * boundvalue) * max(mflux, 0.0)

                        b[id] += b2

                    end


                    if (mflux < 0) && (i != (mesh.l1-1)) && (phi.onoff[i+2,j])
                        num = ((mesh.x[i] + 0.5 * mesh.dx[i]) - mesh.x[i+1])
                        den = (mesh.x[i+1] - mesh.x[i+2])
                        alpha = abs(num / den)
                        aux_mflux = -1.0 * mflux
                        ae, aee, b2 = _convection_secondorderupwind_neighbors_(aux_mflux, alpha)
                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i+1,j]
                        AV[n] = -1.0 * ae

                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i+2,j]
                        AV[n] = -1.0 * aee

                        b[id] += b2

                    elseif (mflux < 0) && (i == (mesh.l1-1))
                        num = ((mesh.x[i] + 0.5 * mesh.dx[i]) - mesh.x[i+1])
                        den = (mesh.x[i+1] - (mesh.x[i+1] + 0.5 * mesh.dx[i+1]))
                        alpha = abs(num / den)
                        aux_mflux = -1.0 * mflux

                        if phi.bounds[i+1,j]
                            boundvalue = find_bondValue((i+1), j, (phi.gIndex[i+1,j]), phi, mesh, bounds, 'e')
                        else
                            boundvalue = phi.eval[i+1,j]
                        end

                        ae = (1.0 + alpha) * max(aux_mflux, 0.0)
                        aee = 0.0
                        b2 = (alpha * boundvalue) * max(aux_mflux, 0.0)

                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i+1,j]
                        AV[n] = -1.0 * ae

                        b[id] += -1.0 * b2

                    elseif (mflux < 0) && (!phi.onoff[i+2,j])
                        num = ((mesh.x[i] + 0.5 * mesh.dx[i]) - mesh.x[i+1])
                        den = (mesh.x[i+1] - (mesh.x[i+1] + 0.5 * mesh.dx[i+1]))
                        alpha = abs(num / den)
                        aux_mflux = -1.0 * mflux

                        if phi.bounds[i+1,j]
                            boundvalue = find_bondValue((i+1), j, (phi.gIndex[i+1,j]), phi, mesh, bounds, 'e')
                        else
                            boundvalue = phi.eval[i+2,j]
                        end

                        ae = (1.0 + alpha) * max(aux_mflux, 0.0)
                        aee = 0.0
                        b2 = (alpha * boundvalue) * max(aux_mflux, 0.0)

                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i+1,j]
                        AV[n] = -1.0 * ae

                        b[id] += -1.0 * b2

                    end

                end

                #South Coefficents
                if (j != 1) && (mesh.m1 != 1) && (phi.onoff[i,j-1])
                    rho = density_interpolation(
                        mesh.dy[j], mesh.dy[j-1], material.ρ[i,j], material.ρ[i,j-1];
                        interpolation = interpolation
                    )

                    mflux = -1.0 * rho * velocityV[i,j] * (mesh.dx[i])
                    mfluxs = -1.0 * rho * velocityV[i,j] * (mesh.dx[i])

                    if (mflux >= 0) && (j != mesh.m1) && (phi.onoff[i,j+1])
                        num = ((mesh.y[j] - 0.5 * mesh.dy[j]) - mesh.y[j])
                        den = (mesh.y[j] - mesh.y[j+1])
                        alpha = abs(num / den)
                        asc, asn, b3 = _convection_secondorderupwind_neighbors_(mflux, alpha)
                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i,j+1]
                        AV[n] = asn
                        b[id] += b3

                    elseif (mflux >= 0) && (j == mesh.m1)
                        num = ((mesh.y[j] - 0.5 * mesh.dy[j]) - mesh.y[j])
                        den = (mesh.y[j] - (mesh.y[j] + 0.5 * mesh.dy[j]))
                        alpha = abs(num / den)

                        if phi.bounds[i,j]
                            boundvalue = find_bondValue(i, j, id, phi, mesh, bounds, 'n')
                        else
                            boundvalue = phi.eval[i,j]
                        end

                        asc = (1.0 + alpha) * max(mflux, 0.0)
                        acfs = (1.0 + alpha) * max(mflux, 0.0)
                        asn = 0.0
                        b3 = (alpha * boundvalue) * max(mflux, 0.0)

                        b[id] += b3

                    elseif (mflux >= 0) && (!phi.onoff[i,j+1])
                        num = ((mesh.y[j] - 0.5 * mesh.dy[j]) - mesh.y[j])
                        den = (mesh.y[j] - (mesh.y[j] + 0.5 * mesh.dy[j]))
                        alpha = abs(num / den)

                        if phi.bounds[i,j]
                            boundvalue = find_bondValue(i, j, id, phi, mesh, bounds, 'n')
                        else
                            boundvalue = phi.eval[i,j+1]
                        end

                        asc = (1.0 + alpha) * max(mflux, 0.0)
                        acfs = (1.0 + alpha) * max(mflux, 0.0)
                        asn = 0.0
                        b3 = (alpha * boundvalue) * max(mflux, 0.0)

                        b[id] += b3

                    end


                    if (mflux < 0) && (j != 2) && (phi.onoff[i,j-2])
                        num = ((mesh.y[j] - 0.5 * mesh.dy[j]) - mesh.y[j-1])
                        den = (mesh.y[j-1] - mesh.y[j-2])
                        alpha = abs(num / den)
                        aux_mflux = -1.0 * mflux
                        as, ass, b3 = _convection_secondorderupwind_neighbors_(aux_mflux, alpha)
                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i,j-1]
                        AV[n] = -1.0 * as

                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i,j-2]
                        AV[n] = -1.0 * ass

                        b[id] += b3

                    elseif (mflux < 0) && (j == 2)
                        num = ((mesh.y[j] - 0.5 * mesh.dy[j]) - mesh.y[j-1])
                        den = (mesh.y[j-1] - (mesh.y[j-1] - 0.5 * mesh.dy[j-1]))
                        alpha = abs(num / den)
                        aux_mflux = -1.0 * mflux

                        if phi.bounds[i,j-1]
                            boundvalue = find_bondValue(i, (j-1), (phi.gIndex[i,j-1]), phi, mesh, bounds, 's')
                        else
                            boundvalue = phi.eval[i,j-1]
                        end

                        as = (1.0 + alpha) * max(aux_mflux, 0.0)
                        ass = 0.0
                        b3 = (alpha * boundvalue) * max(aux_mflux, 0.0)

                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i,j-1]
                        AV[n] = -1.0 * as

                        b[id] += -1.0 * b3

                    elseif (mflux < 0) && (!phi.onoff[i,j-2])
                        num = ((mesh.y[j] - 0.5 * mesh.dy[j]) - mesh.y[j-1])
                        den = (mesh.y[j-1] - (mesh.y[j-1] - 0.5 * mesh.dy[j-1]))
                        alpha = abs(num / den)
                        aux_mflux = -1.0 * mflux

                        if phi.bounds[i,j-1]
                            boundvalue = find_bondValue(i, (j-1), (phi.gIndex[i,j-1]), phi, mesh, bounds, 's')
                        else
                            boundvalue = phi.eval[i,j-2]
                        end

                        as = (1.0 + alpha) * max(aux_mflux, 0.0)
                        ass = 0.0
                        b3 = (alpha * boundvalue) * max(aux_mflux, 0.0)

                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i,j-1]
                        AV[n] = -1.0 * as

                        b[id] += -1.0 * b3

                    end

                end

                #North Coefficents
                if (j != mesh.m1) && (mesh.m1 != 1) && (phi.onoff[i,j+1])
                    rho = density_interpolation(
                        mesh.dy[j], mesh.dy[j+1], material.ρ[i,j], material.ρ[i,j+1];
                        interpolation = interpolation
                    )

                    mflux = rho * velocityV[i,j+1] * (mesh.dx[i])
                    mfluxn = rho * velocityV[i,j+1] * (mesh.dx[i])

                    if (mflux >= 0) && (j != 1) && (phi.onoff[i,j-1])
                        num = ((mesh.y[j] + 0.5 * mesh.dy[j]) - mesh.y[j])
                        den = (mesh.y[j] - mesh.y[j-1])
                        alpha = abs(num / den)
                        anc, ans, b4 = _convection_secondorderupwind_neighbors_(mflux, alpha)
                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i,j-1]
                        AV[n] = ans
                        b[id] += b4

                    elseif (mflux >= 0) && (j == 1)
                        num = ((mesh.y[j] + 0.5 * mesh.dy[j]) - mesh.y[j])
                        den = (mesh.y[j] - (mesh.y[j] - 0.5 * mesh.dy[j]))
                        alpha = abs(num / den)

                        if phi.bounds[i,j]
                            boundvalue = find_bondValue(i, j, id, phi, mesh, bounds, 's')
                        else
                            boundvalue = phi.eval[i,j]
                        end

                        anc = (1.0 + alpha) * max(mflux, 0.0)
                        acfn = (1.0 + alpha) * max(mflux, 0.0)
                        ans = 0.0
                        b4 = (alpha * boundvalue) * max(mflux, 0.0)

                        b[id] += b4

                    elseif (mflux >= 0) && (!phi.onoff[i,j-1])
                        num = ((mesh.y[j] + 0.5 * mesh.dy[j]) - mesh.y[j])
                        den = (mesh.y[j] - (mesh.y[j] - 0.5 * mesh.dy[j]))
                        alpha = abs(num / den)

                        if phi.bounds[i,j]
                            boundvalue = find_bondValue(i, j, id, phi, mesh, bounds, 's')
                        else
                            boundvalue = phi.eval[i,j-1]
                        end

                        anc = (1.0 + alpha) * max(mflux, 0.0)
                        acfn = (1.0 + alpha) * max(mflux, 0.0)
                        ans = 0.0
                        b4 = (alpha * boundvalue) * max(mflux, 0.0)

                        b[id] += b4

                    end


                    if (mflux < 0) && (j != (mesh.m1-1)) && (phi.onoff[i,j+2])
                        num = ((mesh.y[j] + 0.5 * mesh.dy[j]) - mesh.y[j+1])
                        den = (mesh.y[j+1] - mesh.y[j+2])
                        alpha = abs(num / den)
                        aux_mflux = -1.0 * mflux
                        an, ann, b4 = _convection_secondorderupwind_neighbors_(aux_mflux, alpha)
                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i,j+1]
                        AV[n] = -1.0 * an

                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i,j+2]
                        AV[n] = -1.0 * ann

                        b[id] += b4

                    elseif (mflux < 0) && (j == (mesh.m1-1))
                        num = ((mesh.y[j] + 0.5 * mesh.dy[j]) - mesh.y[j+1])
                        den = (mesh.y[j+1] - (mesh.y[j+1] + 0.5 * mesh.dy[j+1]))
                        alpha = abs(num / den)
                        aux_mflux = -1.0 * mflux

                        if phi.bounds[i,j+1]
                            boundvalue = find_bondValue(i, (j+1), (phi.gIndex[i,j+1]), phi, mesh, bounds, 'n')
                        else
                            boundvalue = phi.eval[i,j+1]
                        end

                        an = (1.0 + alpha) * max(aux_mflux, 0.0)
                        ann = 0.0
                        b4 = (alpha * boundvalue) * max(aux_mflux, 0.0)

                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i,j+1]
                        AV[n] = -1.0 * an

                        b[id] += -1.0 * b4

                    elseif (mflux < 0) && (!phi.onoff[i,j+2])
                        num = ((mesh.y[j] + 0.5 * mesh.dy[j]) - mesh.y[j+1])
                        den = (mesh.y[j+1] - (mesh.y[j+1] + 0.5 * mesh.dy[j+1]))
                        alpha = abs(num / den)
                        aux_mflux = -1.0 * mflux

                        if phi.bounds[i,j+1]
                            boundvalue = find_bondValue(i, (j+1), (phi.gIndex[i,j+1]), phi, mesh, bounds, 'n')
                        else
                            boundvalue = phi.eval[i,j+2]
                        end

                        an = (1.0 + alpha) * max(aux_mflux, 0.0)
                        ann = 0.0
                        b4 = (alpha * boundvalue) * max(aux_mflux, 0.0)

                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i,j+1]
                        AV[n] = -1.0 * an

                        b[id] += -1.0 * b4

                    end

                end

                #Center Coefficent
                if phi.bounds[i,j]
                    @inbounds ac, b0 = _convection_secondorderupwind_central_(
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
                    if inout
                        AV[n] = ac + (aec + awc + anc + asc)
                        #AV[n] = ac + (ae + aee + aw + aww + an + ann + as + ass) - (aew + awe + ans + asn) + (acfw + acfe + acfs + acfn) + (mfluxw + mfluxe + mfluxs + mfluxn)
                    else
                        AV[n] = ac + (ae + aee + aw + aww + an + ann + as + ass) - (aew + awe + ans + asn) #ac + (aec + awc + anc + asc)
                    end
                    b[id] += b0
                else
                    n += 1
                    AI[n] = id
                    AJ[n] = id
                    if inout
                        AV[n] = ac + (aec + awc + anc + asc)
                        #AV[n] = ac + (ae + aee + aw + aww + an + ann + as + ass) - (aew + awe + ans + asn) + (acfw + acfe + acfs + acfn) + (mfluxw + mfluxe + mfluxs + mfluxn)
                    else
                        AV[n] = ac + (ae + aee + aw + aww + an + ann + as + ass) - (aew + awe + ans + asn) #ac + (aec + awc + anc + asc)
                    end
                end
            end
        end
    end

    A = sparse(AI[1:n], AJ[1:n], AV[1:n])

    return A, b
end

function _discretize_convection_secondorderupwind_(
    vel::CSVelocity3D,
    phi::CSPhi3D,
    bounds::Dict{String,BoundsStructured},
    material::CSMaterial3D,
    mesh::UnionCSMesh3D,
    inout::Bool = false;
    velocityU::Array{<:AbstractFloat,3} = vel.fValues.uFace,
    velocityV::Array{<:AbstractFloat,3} = vel.fValues.vFace,
    velocityW::Array{<:AbstractFloat,3} = vel.fValues.wFace,
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
    interpolation::Signed = 1,
)
    n_equations = maximum_globalIndex(phi)

    n = 0
    AI = zeros(N, (13 * n_equations))
    AJ = zeros(N, (13 * n_equations))
    AV = zeros(T, (13 * n_equations))

    b = zeros(T, n_equations)

    for i in 1:mesh.l1
        for j in 1:mesh.m1
            for k in 1:mesh.n1
                if phi.onoff[i,j,k]
                    id = phi.gIndex[i,j,k]

                    #Auxiliar variables
                    ac = 0.0
                    b0 = 0.0
                    aw = 0.0
                    awe = 0.0
                    aww = 0.0
                    awc = 0.0
                    b1 = 0.0
                    ae = 0.0
                    aew = 0.0
                    aee = 0.0
                    aec = 0.0
                    b2 = 0.0
                    as = 0.0
                    asn = 0.0
                    ass = 0.0
                    asc = 0.0
                    b3 = 0.0
                    an = 0.0
                    ans = 0.0
                    ann = 0.0
                    anc = 0.0
                    b4 = 0.0
                    ab = 0.0
                    abt = 0.0
                    abb = 0.0
                    abc = 0.0
                    b5 = 0.0
                    at = 0.0
                    atb = 0.0
                    att = 0.0
                    atc = 0.0
                    b6 = 0.0
                    mflux = 0.0
                    num = 0.0
                    den = 0.0

                    #West Coefficents
                    if (i != 1) && (mesh.l1 != 1) && (phi.onoff[i-1,j,k])
                        rho = density_interpolation(
                            mesh.dx[i], mesh.dx[i-1], material.ρ[i,j,k], material.ρ[i-1,j,k];
                            interpolation = interpolation
                        )

                        mflux = -1.0 * rho * velocityU[i,j,k] * (mesh.dy[j] * mesh.dz[k])

                        if (mflux >= 0) && (i != mesh.l1) && (phi.onoff[i+1,j,k])
                            num = ((mesh.x[i] - 0.5 * mesh.dx[i]) - mesh.x[i])
                            den = (mesh.x[i] - mesh.x[i+1])
                            alpha = abs(num / den)
                            awc, awe, b1 = _convection_secondorderupwind_neighbors_(mflux, alpha)
                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i+1,j,k]
                            AV[n] = awe
                            b[id] += b1

                        elseif (mflux >= 0) && (i == mesh.l1)
                            num = ((mesh.x[i] - 0.5 * mesh.dx[i]) - mesh.x[i])
                            den = (mesh.x[i] - (mesh.x[i] + 0.5 * mesh.dx[i]))
                            alpha = abs(num / den)

                            if phi.bounds[i,j,k]
                                boundvalue = find_bondValue(i, j, k, id, phi, mesh, bounds, 'e')
                            else
                                boundvalue = phi.eval[i,j,k]
                            end

                            awc = (1.0 + alpha) * max(mflux, 0.0)
                            awe = 0.0
                            b1 = (alpha * boundvalue) * max(mflux, 0.0)

                            b[id] += b1

                        elseif (mflux >= 0) && (!phi.onoff[i+1,j,k])
                            num = ((mesh.x[i] - 0.5 * mesh.dx[i]) - mesh.x[i])
                            den = (mesh.x[i] - (mesh.x[i] + 0.5 * mesh.dx[i]))
                            alpha = abs(num / den)

                            if phi.bounds[i,j,k]
                                boundvalue = find_bondValue(i, j, k, id, phi, mesh, bounds, 'e')
                            else
                                boundvalue = phi.eval[i+1,j,k]
                            end

                            awc = (1.0 + alpha) * max(mflux, 0.0)
                            awe = 0.0
                            b1 = (alpha * boundvalue) * max(mflux, 0.0)

                            b[id] += b1

                        end


                        if (mflux < 0) && (i != 2) && (phi.onoff[i-2,j,k])
                            num = ((mesh.x[i] - 0.5 * mesh.dx[i]) - mesh.x[i-1])
                            den = (mesh.x[i-1] - mesh.x[i-2])
                            alpha = abs(num / den)
                            aux_mflux = -1.0 * mflux
                            aw, aww, b1 = _convection_secondorderupwind_neighbors_(aux_mflux, alpha)
                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i-1,j,k]
                            AV[n] = -1.0 * aw

                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i-2,j,k]
                            AV[n] = -1.0 * aww

                            b[id] += b1

                        elseif (mflux < 0) && (i == 2)
                            num = ((mesh.x[i] - 0.5 * mesh.dx[i]) - mesh.x[i-1])
                            den = (mesh.x[i-1] - (mesh.x[i-1] - 0.5 * mesh.dx[i-1]))
                            alpha = abs(num / den)
                            aux_mflux = -1.0 * mflux

                            if phi.bounds[i-1,j,k]
                                boundvalue = find_bondValue((i-1), j, k, (phi.gIndex[i-1,j,k]), phi, mesh, bounds, 'w')
                            else
                                boundvalue = phi.eval[i-1,j,k]
                            end

                            aw = (1.0 + alpha) * max(aux_mflux, 0.0)
                            aww = 0.0
                            b1 = (alpha * boundvalue) * max(aux_mflux, 0.0)

                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i-1,j,k]
                            AV[n] = -1.0 * aw

                            b[id] += -1.0 * b1

                        elseif (mflux < 0) && (!phi.onoff[i-2,j,k])
                            num = ((mesh.x[i] - 0.5 * mesh.dx[i]) - mesh.x[i-1])
                            den = (mesh.x[i-1] - (mesh.x[i-1] - 0.5 * mesh.dx[i-1]))
                            alpha = abs(num / den)
                            aux_mflux = -1.0 * mflux

                            if phi.bounds[i-1,j,k]
                                boundvalue = find_bondValue((i-1), j, k, (phi.gIndex[i-1,j,k]), phi, mesh, bounds, 'w')
                            else
                                boundvalue = phi.eval[i-2,j,k]
                            end

                            aw = (1.0 + alpha) * max(aux_mflux, 0.0)
                            aww = 0.0
                            b1 = (alpha * boundvalue) * max(aux_mflux, 0.0)

                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i-1,j,k]
                            AV[n] = -1.0 * aw

                            b[id] += -1.0 * b1

                        end

                    end

                    #East Coefficents
                    if (i != mesh.l1) && (mesh.l1 != 1) && (phi.onoff[i+1,j,k])
                        rho = density_interpolation(
                            mesh.dx[i], mesh.dx[i+1], material.ρ[i,j,k], material.ρ[i+1,j,k];
                            interpolation = interpolation
                        )

                        mflux = rho * velocityU[i+1,j,k] * (mesh.dy[j] * mesh.dz[k])

                        if (mflux >= 0) && (i != 1) && (phi.onoff[i-1,j,k])
                            num = ((mesh.x[i] + 0.5 * mesh.dx[i]) - mesh.x[i])
                            den = (mesh.x[i] - mesh.x[i-1])
                            alpha = abs(num / den)
                            aec, aew, b2 = _convection_secondorderupwind_neighbors_(mflux, alpha)
                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i-1,j,k]
                            AV[n] = aew
                            b[id] += b2

                        elseif (mflux >= 0) && (i == 1)
                            num = ((mesh.x[i] + 0.5 * mesh.dx[i]) - mesh.x[i])
                            den = (mesh.x[i] - (mesh.x[i] - 0.5 * mesh.dx[i]))
                            alpha = abs(num / den)

                            if phi.bounds[i,j,k]
                                boundvalue = find_bondValue(i, j, k, id, phi, mesh, bounds, 'w')
                            else
                                boundvalue = phi.eval[i,j,k]
                            end

                            aec = (1.0 + alpha) * max(mflux, 0.0)
                            aew = 0.0
                            b2 = (alpha * boundvalue) * max(mflux, 0.0)

                            b[id] += b2

                        elseif (mflux >= 0) && (!phi.onoff[i-1,j,k])
                            num = ((mesh.x[i] + 0.5 * mesh.dx[i]) - mesh.x[i])
                            den = (mesh.x[i] - (mesh.x[i] - 0.5 * mesh.dx[i]))
                            alpha = abs(num / den)

                            if phi.bounds[i,j,k]
                                boundvalue = find_bondValue(i, j, k, id, phi, mesh, bounds, 'w')
                            else
                                boundvalue = phi.eval[i-1,j,k]
                            end

                            aec = (1.0 + alpha) * max(mflux, 0.0)
                            aew = 0.0
                            b2 = (alpha * boundvalue) * max(mflux, 0.0)

                            b[id] += b2

                        end


                        if (mflux < 0) && (i != (mesh.l1-1)) && (phi.onoff[i+2,j,k])
                            num = ((mesh.x[i] + 0.5 * mesh.dx[i]) - mesh.x[i+1])
                            den = (mesh.x[i+1] - mesh.x[i+2])
                            alpha = abs(num / den)
                            aux_mflux = -1.0 * mflux
                            ae, aee, b2 = _convection_secondorderupwind_neighbors_(aux_mflux, alpha)
                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i+1,j,k]
                            AV[n] = -1.0 * ae

                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i+2,j,k]
                            AV[n] = -1.0 * aee

                            b[id] += b2

                        elseif (mflux < 0) && (i == (mesh.l1-1))
                            num = ((mesh.x[i] + 0.5 * mesh.dx[i]) - mesh.x[i+1])
                            den = (mesh.x[i+1] - (mesh.x[i+1] + 0.5 * mesh.dx[i+1]))
                            alpha = abs(num / den)
                            aux_mflux = -1.0 * mflux

                            if phi.bounds[i+1,j,k]
                                boundvalue = find_bondValue((i+1), j, k, (phi.gIndex[i+1,j,k]), phi, mesh, bounds, 'e')
                            else
                                boundvalue = phi.eval[i+1,j,k]
                            end

                            ae = (1.0 + alpha) * max(aux_mflux, 0.0)
                            aee = 0.0
                            b2 = (alpha * boundvalue) * max(aux_mflux, 0.0)

                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i+1,j,k]
                            AV[n] = -1.0 * ae

                            b[id] += -1.0 * b2

                        elseif (mflux < 0) && (!phi.onoff[i+2,j,k])
                            num = ((mesh.x[i] + 0.5 * mesh.dx[i]) - mesh.x[i+1])
                            den = (mesh.x[i+1] - (mesh.x[i+1] + 0.5 * mesh.dx[i+1]))
                            alpha = abs(num / den)
                            aux_mflux = -1.0 * mflux

                            if phi.bounds[i+1,j,k]
                                boundvalue = find_bondValue((i+1), j, k, (phi.gIndex[i+1,j,k]), phi, mesh, bounds, 'e')
                            else
                                boundvalue = phi.eval[i+2,j,k]
                            end

                            ae = (1.0 + alpha) * max(aux_mflux, 0.0)
                            aee = 0.0
                            b2 = (alpha * boundvalue) * max(aux_mflux, 0.0)

                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i+1,j,k]
                            AV[n] = -1.0 * ae

                            b[id] += -1.0 * b2

                        end

                    end

                    #South Coefficents
                    if (j != 1) && (mesh.m1 != 1) && (phi.onoff[i,j-1,k])
                        rho = density_interpolation(
                            mesh.dy[j], mesh.dy[j-1], material.ρ[i,j,k], material.ρ[i,j-1,k];
                            interpolation = interpolation
                        )

                        mflux = -1.0 * rho * velocityV[i,j,k] * (mesh.dx[i] * mesh.dz[k])

                        if (mflux >= 0) && (j != mesh.m1) && (phi.onoff[i,j+1,k])
                            num = ((mesh.y[j] - 0.5 * mesh.dy[j]) - mesh.y[j])
                            den = (mesh.y[j] - mesh.y[j+1])
                            alpha = abs(num / den)
                            asc, asn, b3 = _convection_secondorderupwind_neighbors_(mflux, alpha)
                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i,j+1,k]
                            AV[n] = asn
                            b[id] += b3

                        elseif (mflux >= 0) && (j == mesh.m1)
                            num = ((mesh.y[j] - 0.5 * mesh.dy[j]) - mesh.y[j])
                            den = (mesh.y[j] - (mesh.y[j] + 0.5 * mesh.dy[j]))
                            alpha = abs(num / den)

                            if phi.bounds[i,j,k]
                                boundvalue = find_bondValue(i, j, k, id, phi, mesh, bounds, 'n')
                            else
                                boundvalue = phi.eval[i,j,k]
                            end

                            asc = (1.0 + alpha) * max(mflux, 0.0)
                            asn = 0.0
                            b3 = (alpha * boundvalue) * max(mflux, 0.0)

                            b[id] += b3

                        elseif (mflux >= 0) && (!phi.onoff[i,j+1,k])
                            num = ((mesh.y[j] - 0.5 * mesh.dy[j]) - mesh.y[j])
                            den = (mesh.y[j] - (mesh.y[j] + 0.5 * mesh.dy[j]))
                            alpha = abs(num / den)

                            if phi.bounds[i,j,k]
                                boundvalue = find_bondValue(i, j, k, id, phi, mesh, bounds, 'n')
                            else
                                boundvalue = phi.eval[i,j+1,k]
                            end

                            asc = (1.0 + alpha) * max(mflux, 0.0)
                            asn = 0.0
                            b3 = (alpha * boundvalue) * max(mflux, 0.0)

                            b[id] += b3

                        end


                        if (mflux < 0) && (j != 2) && (phi.onoff[i,j-2,k])
                            num = ((mesh.y[j] - 0.5 * mesh.dy[j]) - mesh.y[j-1])
                            den = (mesh.y[j-1] - mesh.y[j-2])
                            alpha = abs(num / den)
                            aux_mflux = -1.0 * mflux
                            as, ass, b3 = _convection_secondorderupwind_neighbors_(aux_mflux, alpha)
                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i,j-1,k]
                            AV[n] = -1.0 * as

                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i,j-2,k]
                            AV[n] = -1.0 * ass

                            b[id] += b3

                        elseif (mflux < 0) && (j == 2)
                            num = ((mesh.y[j] - 0.5 * mesh.dy[j]) - mesh.y[j-1])
                            den = (mesh.y[j-1] - (mesh.y[j-1] - 0.5 * mesh.dy[j-1]))
                            alpha = abs(num / den)
                            aux_mflux = -1.0 * mflux

                            if phi.bounds[i,j-1,k]
                                boundvalue = find_bondValue(i, (j-1), k, (phi.gIndex[i,j-1,k]), phi, mesh, bounds, 's')
                            else
                                boundvalue = phi.eval[i,j-1,k]
                            end

                            as = (1.0 + alpha) * max(aux_mflux, 0.0)
                            ass = 0.0
                            b3 = (alpha * boundvalue) * max(aux_mflux, 0.0)

                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i,j-1,k]
                            AV[n] = -1.0 * as

                            b[id] += -1.0 * b3

                        elseif (mflux < 0) && (!phi.onoff[i,j-2,k])
                            num = ((mesh.y[j] - 0.5 * mesh.dy[j]) - mesh.y[j-1])
                            den = (mesh.y[j-1] - (mesh.y[j-1] - 0.5 * mesh.dy[j-1]))
                            alpha = abs(num / den)
                            aux_mflux = -1.0 * mflux

                            if phi.bounds[i,j-1,k]
                                boundvalue = find_bondValue(i, (j-1), k, (phi.gIndex[i,j-1,k]), phi, mesh, bounds, 's')
                            else
                                boundvalue = phi.eval[i,j-2,k]
                            end

                            as = (1.0 + alpha) * max(aux_mflux, 0.0)
                            ass = 0.0
                            b3 = (alpha * boundvalue) * max(aux_mflux, 0.0)

                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i,j-1,k]
                            AV[n] = -1.0 * as

                            b[id] += -1.0 * b3

                        end

                    end

                    #North Coefficents
                    if (j != mesh.m1) && (mesh.m1 != 1) && (phi.onoff[i,j+1,k])
                        rho = density_interpolation(
                            mesh.dy[j], mesh.dy[j+1], material.ρ[i,j,k], material.ρ[i,j+1,k];
                            interpolation = interpolation
                        )

                        mflux = rho * velocityV[i,j+1,k] * (mesh.dx[i] * mesh.dz[k])

                        if (mflux >= 0) && (j != 1) && (phi.onoff[i,j-1,k])
                            num = ((mesh.y[j] + 0.5 * mesh.dy[j]) - mesh.y[j])
                            den = (mesh.y[j] - mesh.y[j-1])
                            alpha = abs(num / den)
                            anc, ans, b4 = _convection_secondorderupwind_neighbors_(mflux, alpha)
                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i,j-1,k]
                            AV[n] = ans
                            b[id] += b4

                        elseif (mflux >= 0) && (j == 1)
                            num = ((mesh.y[j] + 0.5 * mesh.dy[j]) - mesh.y[j])
                            den = (mesh.y[j] - (mesh.y[j] - 0.5 * mesh.dy[j]))
                            alpha = abs(num / den)

                            if phi.bounds[i,j,k]
                                boundvalue = find_bondValue(i, j, k, id, phi, mesh, bounds, 's')
                            else
                                boundvalue = phi.eval[i,j,k]
                            end

                            anc = (1.0 + alpha) * max(mflux, 0.0)
                            ans = 0.0
                            b4 = (alpha * boundvalue) * max(mflux, 0.0)

                            b[id] += b4

                        elseif (mflux >= 0) && (!phi.onoff[i,j-1,k])
                            num = ((mesh.y[j] + 0.5 * mesh.dy[j]) - mesh.y[j])
                            den = (mesh.y[j] - (mesh.y[j] - 0.5 * mesh.dy[j]))
                            alpha = abs(num / den)

                            if phi.bounds[i,j,k]
                                boundvalue = find_bondValue(i, j, k, id, phi, mesh, bounds, 's')
                            else
                                boundvalue = phi.eval[i,j-1,k]
                            end

                            anc = (1.0 + alpha) * max(mflux, 0.0)
                            ans = 0.0
                            b4 = (alpha * boundvalue) * max(mflux, 0.0)

                            b[id] += b4

                        end


                        if (mflux < 0) && (j != (mesh.m1-1)) && (phi.onoff[i,j+2,k])
                            num = ((mesh.y[j] + 0.5 * mesh.dy[j]) - mesh.y[j+1])
                            den = (mesh.y[j+1] - mesh.y[j+2])
                            alpha = abs(num / den)
                            aux_mflux = -1.0 * mflux
                            an, ann, b4 = _convection_secondorderupwind_neighbors_(aux_mflux, alpha)
                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i,j+1,k]
                            AV[n] = -1.0 * an

                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i,j+2,k]
                            AV[n] = -1.0 * ann

                            b[id] += b4

                        elseif (mflux < 0) && (j == (mesh.m1-1))
                            num = ((mesh.y[j] + 0.5 * mesh.dy[j]) - mesh.y[j+1])
                            den = (mesh.y[j+1] - (mesh.y[j+1] + 0.5 * mesh.dy[j+1]))
                            alpha = abs(num / den)
                            aux_mflux = -1.0 * mflux

                            if phi.bounds[i,j+1,k]
                                boundvalue = find_bondValue(i, (j+1), k, (phi.gIndex[i,j+1,k]), phi, mesh, bounds, 'n')
                            else
                                boundvalue = phi.eval[i,j+1,k]
                            end

                            an = (1.0 + alpha) * max(aux_mflux, 0.0)
                            ann = 0.0
                            b4 = (alpha * boundvalue) * max(aux_mflux, 0.0)

                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i,j+1,k]
                            AV[n] = -1.0 * an

                            b[id] += -1.0 * b4

                        elseif (mflux < 0) && (!phi.onoff[i,j+2,k])
                            num = ((mesh.y[j] + 0.5 * mesh.dy[j]) - mesh.y[j+1])
                            den = (mesh.y[j+1] - (mesh.y[j+1] + 0.5 * mesh.dy[j+1]))
                            alpha = abs(num / den)
                            aux_mflux = -1.0 * mflux

                            if phi.bounds[i,j+1,k]
                                boundvalue = find_bondValue(i, (j+1), k, (phi.gIndex[i,j+1,k]), phi, mesh, bounds, 'n')
                            else
                                boundvalue = phi.eval[i,j+2,k]
                            end

                            an = (1.0 + alpha) * max(aux_mflux, 0.0)
                            ann = 0.0
                            b4 = (alpha * boundvalue) * max(aux_mflux, 0.0)

                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i,j+1,k]
                            AV[n] = -1.0 * an

                            b[id] += -1.0 * b4

                        end

                    end

                    #Bottom Coefficents
                    if (k != 1) && (mesh.n1 != 1) && (phi.onoff[i,j,k-1])
                        rho = density_interpolation(
                            mesh.dz[k], mesh.dz[k-1], material.ρ[i,j,k], material.ρ[i,j,k-1];
                            interpolation = interpolation
                        )

                        mflux = -1.0 * rho * velocityW[i,j,k] * (mesh.dx[i] * mesh.dy[j])

                        if (mflux >= 0) && (k != mesh.n1) && (phi.onoff[i,j,k+1])
                            num = ((mesh.z[k] - 0.5 * mesh.dz[k]) - mesh.z[k])
                            den = (mesh.z[k] - mesh.z[k+1])
                            alpha = abs(num / den)
                            abc, abt, b5 = _convection_secondorderupwind_neighbors_(mflux, alpha)
                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i,j,k+1]
                            AV[n] = abt
                            b[id] += b5

                        elseif (mflux >= 0) && (k == mesh.n1)
                            num = ((mesh.z[k] - 0.5 * mesh.dz[k]) - mesh.z[k])
                            den = (mesh.z[k] - (mesh.z[k] + 0.5 * mesh.dz[k]))
                            alpha = abs(num / den)

                            if phi.bounds[i,j,k]
                                boundvalue = find_bondValue(i, j, k, id, phi, mesh, bounds, 't')
                            else
                                boundvalue = phi.eval[i,j,k]
                            end

                            abc = (1.0 + alpha) * max(mflux, 0.0)
                            abt = 0.0
                            b5 = (alpha * boundvalue) * max(mflux, 0.0)

                            b[id] += b5

                        elseif (mflux >= 0) && (!phi.onoff[i,j,k+1])
                            num = ((mesh.z[k] - 0.5 * mesh.dz[k]) - mesh.z[k])
                            den = (mesh.z[k] - (mesh.z[k] + 0.5 * mesh.dz[k]))
                            alpha = abs(num / den)

                            if phi.bounds[i,j,k]
                                boundvalue = find_bondValue(i, j, k, id, phi, mesh, bounds, 't')
                            else
                                boundvalue = phi.eval[i+1,j,k]
                            end

                            abc = (1.0 + alpha) * max(mflux, 0.0)
                            abt = 0.0
                            b5 = (alpha * boundvalue) * max(mflux, 0.0)

                            b[id] += b5

                        end


                        if (mflux < 0) && (k != 2) && (phi.onoff[i,j,k-2])
                            num = ((mesh.z[k] - 0.5 * mesh.dz[k]) - mesh.z[k-1])
                            den = (mesh.z[k-1] - mesh.z[k-2])
                            alpha = abs(num / den)
                            aux_mflux = -1.0 * mflux
                            ab, abb, b5 = _convection_secondorderupwind_neighbors_(aux_mflux, alpha)
                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i,j,k-1]
                            AV[n] = -1.0 * ab

                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i,j,k-2]
                            AV[n] = -1.0 * abb

                            b[id] += b5

                        elseif (mflux < 0) && (k == 2)
                            num = ((mesh.z[k] - 0.5 * mesh.dz[k]) - mesh.z[k-1])
                            den = (mesh.z[k-1] - (mesh.z[k-1] - 0.5 * mesh.dz[k-1]))
                            alpha = abs(num / den)
                            aux_mflux = -1.0 * mflux

                            if phi.bounds[i,j,k-1]
                                boundvalue = find_bondValue(i, j, (k-1), (phi.gIndex[i,j,k-1]), phi, mesh, bounds, 'b')
                            else
                                boundvalue = phi.eval[i,j,k-1]
                            end

                            ab = (1.0 + alpha) * max(aux_mflux, 0.0)
                            abb = 0.0
                            b5 = (alpha * boundvalue) * max(aux_mflux, 0.0)

                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i,j,k-1]
                            AV[n] = -1.0 * ab

                            b[id] += -1.0 * b5

                        elseif (mflux < 0) && (!phi.onoff[i,j,k-2])
                            num = ((mesh.z[k] - 0.5 * mesh.dz[k]) - mesh.z[k-1])
                            den = (mesh.z[k-1] - (mesh.z[k-1] - 0.5 * mesh.dz[k-1]))
                            alpha = abs(num / den)
                            aux_mflux = -1.0 * mflux

                            if phi.bounds[i,j,k-1]
                                boundvalue = find_bondValue(i, j, (k-1), (phi.gIndex[i,j,k-1]), phi, mesh, bounds, 'b')
                            else
                                boundvalue = phi.eval[i,j,k-2]
                            end

                            ab = (1.0 + alpha) * max(aux_mflux, 0.0)
                            abb = 0.0
                            b5 = (alpha * boundvalue) * max(aux_mflux, 0.0)

                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i,j,k-1]
                            AV[n] = -1.0 * ab

                            b[id] += -1.0 * b5

                        end

                    end

                    #Top Coefficents
                    if (k != mesh.n1) && (mesh.n1 != 1) && (phi.onoff[i,j,k+1])
                        rho = density_interpolation(
                            mesh.dz[k], mesh.dz[k+1], material.ρ[i,j,k], material.ρ[i,j,k+1];
                            interpolation = interpolation
                        )

                        mflux = rho * velocityW[i,j,k+1] * (mesh.dx[i] * mesh.dy[j])

                        if (mflux >= 0) && (k != 1) && (phi.onoff[i,j,k-1])
                            num = ((mesh.z[k] + 0.5 * mesh.dz[k]) - mesh.z[k])
                            den = (mesh.z[k] - mesh.z[k-1])
                            alpha = abs(num / den)
                            atc, atb, b6 = _convection_secondorderupwind_neighbors_(mflux, alpha)
                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i,j,k-1]
                            AV[n] = atb
                            b[id] += b6

                        elseif (mflux >= 0) && (k == 1)
                            num = ((mesh.z[k] + 0.5 * mesh.dz[k]) - mesh.z[k])
                            den = (mesh.z[k] - (mesh.z[k] - 0.5 * mesh.dz[k]))
                            alpha = abs(num / den)

                            if phi.bounds[i,j,k]
                                boundvalue = find_bondValue(i, j, k, id, phi, mesh, bounds, 'b')
                            else
                                boundvalue = phi.eval[i,j,k]
                            end

                            atc = (1.0 + alpha) * max(mflux, 0.0)
                            atb = 0.0
                            b6 = (alpha * boundvalue) * max(mflux, 0.0)

                            b[id] += b6

                        elseif (mflux >= 0) && (!phi.onoff[i,j,k-1])
                            num = ((mesh.z[k] + 0.5 * mesh.dz[k]) - mesh.z[k])
                            den = (mesh.z[k] - (mesh.z[k] - 0.5 * mesh.dz[k]))
                            alpha = abs(num / den)

                            if phi.bounds[i,j,k]
                                boundvalue = find_bondValue(i, j, k, id, phi, mesh, bounds, 'b')
                            else
                                boundvalue = phi.eval[i,j,k-1]
                            end

                            atc = (1.0 + alpha) * max(mflux, 0.0)
                            atb = 0.0
                            b6 = (alpha * boundvalue) * max(mflux, 0.0)

                            b[id] += b6

                        end


                        if (mflux < 0) && (k != (mesh.n1-1)) && (phi.onoff[i,j,k+2])
                            num = ((mesh.z[k] + 0.5 * mesh.dz[k]) - mesh.z[k+1])
                            den = (mesh.z[k+1] - mesh.z[k+2])
                            alpha = abs(num / den)
                            aux_mflux = -1.0 * mflux
                            at, att, b6 = _convection_secondorderupwind_neighbors_(aux_mflux, alpha)
                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i,j,k+1]
                            AV[n] = -1.0 * at

                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i,j,k+2]
                            AV[n] = -1.0 * att

                            b[id] += b6

                        elseif (mflux < 0) && (k == (mesh.n1-1))
                            num = ((mesh.z[k] + 0.5 * mesh.dz[k]) - mesh.z[k+1])
                            den = (mesh.z[k+1] - (mesh.z[k+1] + 0.5 * mesh.dz[k+1]))
                            alpha = abs(num / den)
                            aux_mflux = -1.0 * mflux

                            if phi.bounds[i,j,k+1]
                                boundvalue = find_bondValue(i, j, (k+1), (phi.gIndex[i,j,k+1]), phi, mesh, bounds, 't')
                            else
                                boundvalue = phi.eval[i,j,k+1]
                            end

                            at = (1.0 + alpha) * max(aux_mflux, 0.0)
                            att = 0.0
                            b6 = (alpha * boundvalue) * max(aux_mflux, 0.0)

                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i,j,k+1]
                            AV[n] = -1.0 * at

                            b[id] += -1.0 * b6

                        elseif (mflux < 0) && (!phi.onoff[i,j,k+2])
                            num = ((mesh.z[k] + 0.5 * mesh.dz[k]) - mesh.z[k+1])
                            den = (mesh.z[k+1] - (mesh.z[k+1] + 0.5 * mesh.dz[k+1]))
                            alpha = abs(num / den)
                            aux_mflux = -1.0 * mflux

                            if phi.bounds[i,j,k+1]
                                boundvalue = find_bondValue(i, j, (k+1), (phi.gIndex[i,j,k+1]), phi, mesh, bounds, 't')
                            else
                                boundvalue = phi.eval[i,j,k+2]
                            end

                            at = (1.0 + alpha) * max(aux_mflux, 0.0)
                            att = 0.0
                            b6 = (alpha * boundvalue) * max(aux_mflux, 0.0)

                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i,j,k+1]
                            AV[n] = -1.0 * at

                            b[id] += -1.0 * b6

                        end

                    end

                    #Center Coefficent
                    if phi.bounds[i,j,k]
                        @inbounds ac, b0 = _convection_secondorderupwind_central_(
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
                        if inout
                            AV[n] = ac + aec + awc + anc + asc + atc + abc
                        else
                            AV[n] = ac + (ae + aee + aw + aww + an + ann + as + ass + at + att + ab + abb) - (aew + awe + ans + asn + atb + abt) # ac + aec + awc + anc + asc + atc + abc
                        end
                        b[id] += b0
                    else
                        n += 1
                        AI[n] = id
                        AJ[n] = id
                        if inout
                            AV[n] = ac + aec + awc + anc + asc + atc + abc
                        else
                            AV[n] = ac + (ae + aee + aw + aww + an + ann + as + ass + at + att + ab + abb) - (aew + awe + ans + asn + atb + abt) # ac + aec + awc + anc + asc + atc + abc
                        end
                    end
                end
            end
        end
    end

    A = sparse(AI[1:n], AJ[1:n], AV[1:n])

    return A, b
end
