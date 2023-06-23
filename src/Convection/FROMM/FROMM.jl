"""

"""
function _discretize_convection_FROMM_(
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
    # println("TVD")

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
                awc = 0.0
                awd1 = 0.0
                awd2 = 0.0
                aec = 0.0
                aed1 = 0.0
                aed2 = 0.0
                asc = 0.0
                asd1 = 0.0
                asd2 = 0.0
                anc = 0.0
                and1 = 0.0
                and2 = 0.0
                mfluxw = 0.0
                mfluxe = 0.0
                mfluxs = 0.0
                mfluxn = 0.0

                ge = 0.0
                gw = 0.0
                gee = 0.0
                gww = 0.0
                gn = 0.0
                gs = 0.0
                gnn = 0.0
                gss = 0.0

                vu = 0.0
                vd = 0.0
                vc = 0.0

                cf_vu = 0.0
                cf_vd = 0.0
                cf_vc = 0.0

                #West Coefficents
                if (i != 1) && (mesh.l1 != 1) && (phi.onoff[i-1,j])
                    rho = density_interpolation(
                        mesh.dx[i], mesh.dx[i-1], material.ρ[i,j], material.ρ[i-1,j];
                        interpolation = interpolation
                    )

                    mfluxw = -1.0 * rho * velocityU[i,j] * (mesh.dy[j])

                    if (mfluxw >= 0)
                        # Geometric Parameters
                        # ge
                        if (i == mesh.l1) || (!phi.onoff[i+1,j])
                            ge = 1.0
                        else
                            ge = (mesh.dx[i]) / (mesh.dx[i] + mesh.dx[i+1])
                        end
                        # gw
                        gw = (mesh.dx[i]) / (mesh.dx[i] + mesh.dx[i-1])

                        # UPWIND
                        if (i == mesh.l1) || (!phi.onoff[i+1,j])
                            if phi.bounds[i,j]
                                vu = find_bondValue(i, j, id, phi, mesh, bounds, 'e')
                            else
                                vu = phi.eval[i,j]
                            end
                        else
                            vu = phi.eval[i+1,j]
                        end

                        # DOWNWIND
                        vd = phi.eval[i-1,j]

                        # CENTER
                        vc = phi.eval[i,j]

                        m_coef = 1.0
                        n_coef = 1.0
                        c_coef = 0.0

                        m_tvd = m_coef * ge
                        n_tvd = n_coef * gw + c_coef

                        cf_vc = (1.0 + 0.5 * m_tvd - 0.5 * n_tvd) * mfluxw
                        cf_vd = (0.5 * n_tvd) * mfluxw
                        cf_vu = (-0.5 * m_tvd) * mfluxw

                        awc = cf_vc
                        awd1 = cf_vd
                        awd2 = cf_vu

                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i-1,j]
                        AV[n] = awd1

                        if (i == mesh.l1) || (!phi.onoff[i+1,j])
                            b[id] += -1.0 * awd2 * vu
                        else
                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i+1,j]
                            AV[n] = awd2
                        end

                    elseif (mfluxw < 0)
                        # Geometric Parameters
                        # gww
                        if (i == 2) || (!phi.onoff[i-2,j])
                            gww = 1.0
                        else
                            gww = (mesh.dx[i-1]) / (mesh.dx[i-1] + mesh.dx[i-2])
                        end
                        # gw
                        gw = (mesh.dx[i]) / (mesh.dx[i] + mesh.dx[i-1])

                        # UPWIND
                        if (i == 2) || (!phi.onoff[i-2,j])
                            if phi.bounds[i-1,j]
                                vu = find_bondValue(i-1, j, phi.gIndex[i-1,j], phi, mesh, bounds, 'w')
                            else
                                vu = phi.eval[i-1,j]
                            end
                        else
                            vu = phi.eval[i-2,j]
                        end

                        # DOWNWIND
                        vd = phi.eval[i,j]

                        # CENTER
                        vc = phi.eval[i-1,j]

                        m_coef = 1.0
                        n_coef = 1.0
                        c_coef = 0.0

                        m_tvd = m_coef * gww
                        n_tvd = n_coef * (1.0 - gw) + c_coef

                        cf_vc = (1.0 + 0.5 * m_tvd - 0.5 * n_tvd) * mfluxw
                        cf_vd = (0.5 * n_tvd) * mfluxw
                        cf_vu = (-0.5 * m_tvd) * mfluxw

                        awc = cf_vd
                        awd1 = cf_vc
                        awd2 = cf_vu

                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i-1,j]
                        AV[n] = awd1

                        if (i == 2) || (!phi.onoff[i-2,j])
                            b[id] += -1.0 * awd2 * vu
                        else
                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i-2,j]
                            AV[n] = awd2
                        end
                    end
                end

                #East Coefficents
                if (i != mesh.l1) && (mesh.l1 != 1) && (phi.onoff[i+1,j])
                    rho = density_interpolation(
                        mesh.dx[i], mesh.dx[i+1], material.ρ[i,j], material.ρ[i+1,j];
                        interpolation = interpolation
                    )

                    mfluxe = rho * velocityU[i+1,j] * (mesh.dy[j])

                    if (mfluxe >= 0)
                        # Geometric Parameters
                        # gw
                        if (i == 1) || (!phi.onoff[i-1,j])
                            gw = 1.0
                        else
                            gw = (mesh.dx[i]) / (mesh.dx[i] + mesh.dx[i-1])
                        end
                        # ge
                        ge = (mesh.dx[i]) / (mesh.dx[i] + mesh.dx[i+1])

                        # UPWIND
                        if (i == 1) || (!phi.onoff[i-1,j])
                            if phi.bounds[i,j]
                                vu = find_bondValue(i, j, id, phi, mesh, bounds, 'w')
                            else
                                vu = phi.eval[i,j]
                            end
                        else
                            vu = phi.eval[i-1,j]
                        end

                        # DOWNWIND
                        vd = phi.eval[i+1,j]

                        # CENTER
                        vc = phi.eval[i,j]

                        m_coef = 1.0
                        n_coef = 1.0
                        c_coef = 0.0

                        m_tvd = m_coef * gw
                        n_tvd = n_coef * ge + c_coef

                        cf_vc = (1.0 + 0.5 * m_tvd - 0.5 * n_tvd) * mfluxe
                        cf_vd = (0.5 * n_tvd) * mfluxe
                        cf_vu = (-0.5 * m_tvd) * mfluxe

                        aec = cf_vc
                        aed1 = cf_vd
                        aed2 = cf_vu

                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i+1,j]
                        AV[n] = aed1

                        if (i == 1) || (!phi.onoff[i-1,j])
                            b[id] += -1.0 * aed2 * vu
                        else
                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i-1,j]
                            AV[n] = aed2
                        end

                    elseif (mfluxe < 0)
                        # Geometric Parameters
                        # ge
                        ge = (mesh.dx[i]) / (mesh.dx[i] + mesh.dx[i+1])
                        # gee
                        if (i == mesh.l1-1) || (!phi.onoff[i+2,j])
                            gee = 1.0
                        else
                            gee = (mesh.dx[i+1]) / (mesh.dx[i+1] + mesh.dx[i+2])
                        end

                        # UPWIND
                        if (i == mesh.l1-1) || (!phi.onoff[i+2,j])
                            if phi.bounds[i,j]
                                vu = find_bondValue(i+1, j, phi.gIndex[i+1,j], phi, mesh, bounds, 'e')
                            else
                                vu = phi.eval[i+1,j]
                            end
                        else
                            vu = phi.eval[i+2,j]
                        end

                        # DOWNWIND
                        vd = phi.eval[i,j]

                        # CENTER
                        vc = phi.eval[i+1,j]

                        m_coef = 1.0
                        n_coef = 1.0
                        c_coef = 0.0

                        m_tvd = m_coef * gee
                        n_tvd = n_coef * (1.0 - ge) + c_coef

                        cf_vc = (1.0 + 0.5 * m_tvd - 0.5 * n_tvd) * mfluxe
                        cf_vd = (0.5 * n_tvd) * mfluxe
                        cf_vu = (-0.5 * m_tvd) * mfluxe

                        aec = cf_vd
                        aed1 = cf_vc
                        aed2 = cf_vu

                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i+1,j]
                        AV[n] = aed1

                        if (i == mesh.l1-1) || (!phi.onoff[i+2,j])
                            b[id] += -1.0 * aed2 * vu
                        else
                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i+2,j]
                            AV[n] = aed2
                        end

                    end

                end

                # South Coefficents
                if (j != 1) && (mesh.m1 != 1) && (phi.onoff[i,j-1])
                    rho = density_interpolation(
                        mesh.dy[j], mesh.dy[j-1], material.ρ[i,j], material.ρ[i,j-1];
                        interpolation = interpolation
                    )

                    mfluxs = -1.0 * rho * velocityV[i,j] * (mesh.dx[i])

                    if (mfluxs >= 0)
                        # Geometric Parameters
                        # gn
                        if (j == mesh.m1) || (!phi.onoff[i,j+1])
                            gn = 1.0
                        else
                            gn = (mesh.dy[j]) / (mesh.dy[j] + mesh.dy[j+1])
                        end
                        # gs
                        gs = (mesh.dy[j]) / (mesh.dy[j] + mesh.dy[j-1])

                        # UPWIND
                        if (j == mesh.m1) || (!phi.onoff[i,j+1])
                            if phi.bounds[i,j]
                                vu = find_bondValue(i, j, id, phi, mesh, bounds, 'n')
                            else
                                vu = phi.eval[i,j]
                            end
                        else
                            vu = phi.eval[i,j+1]
                        end

                        # DOWNWIND
                        vd = phi.eval[i,j-1]

                        # CENTER
                        vc = phi.eval[i,j]

                        m_coef = 1.0
                        n_coef = 1.0
                        c_coef = 0.0

                        m_tvd = m_coef * gn
                        n_tvd = n_coef * gs + c_coef

                        cf_vc = (1.0 + 0.5 * m_tvd - 0.5 * n_tvd) * mfluxs
                        cf_vd = (0.5 * n_tvd) * mfluxs
                        cf_vu = (-0.5 * m_tvd) * mfluxs

                        asc = cf_vc
                        asd1 = cf_vd
                        asd2 = cf_vu

                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i,j-1]
                        AV[n] = asd1

                        if (j == mesh.m1)
                            b[id] += -1.0 * asd2 * vu
                        else
                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i,j+1]
                            AV[n] = asd2
                        end

                    elseif (mfluxs < 0)
                        # Geometric Parameters
                        # gss
                        if (j == 2) || (!phi.onoff[i,j-2])
                            gss = 1.0
                        else
                            gss = (mesh.dy[j-1]) / (mesh.dy[j-1] + mesh.dy[j-2])
                        end
                        # gs
                        gs = (mesh.dy[j]) / (mesh.dy[j] + mesh.dy[j-1])

                        # UPWIND
                        if (j == 2) || (!phi.onoff[i,j-2])
                            if phi.bounds[i,j-1]
                                vu = find_bondValue(i, j-1, phi.gIndex[i,j-1], phi, mesh, bounds, 's')
                            else
                                vu = phi.eval[i,j-1]
                            end
                        else
                            vu = phi.eval[i,j-2]
                        end

                        # DOWNWIND
                        vd = phi.eval[i,j]

                        # CENTER
                        vc = phi.eval[i,j-1]

                        m_coef = 1.0
                        n_coef = 1.0
                        c_coef = 0.0

                        m_tvd = m_coef * gss
                        n_tvd = n_coef * (1.0 - gs) + c_coef

                        cf_vc = (1.0 + 0.5 * m_tvd - 0.5 * n_tvd) * mfluxs
                        cf_vd = (0.5 * n_tvd) * mfluxs
                        cf_vu = (-0.5 * m_tvd) * mfluxs

                        asc = cf_vd
                        asd1 = cf_vc
                        asd2 = cf_vu

                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i,j-1]
                        AV[n] = asd1

                        if (j == 2) || (!phi.onoff[i,j-2])
                            b[id] += -1.0 * asd2 * vu
                        else
                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i,j-2]
                            AV[n] = asd2
                        end

                    end

                end

                #North Coefficents
                if (j != mesh.m1) && (mesh.m1 != 1) && (phi.onoff[i,j+1])
                    rho = density_interpolation(
                        mesh.dy[j], mesh.dy[j+1], material.ρ[i,j], material.ρ[i,j+1];
                        interpolation = interpolation
                    )

                    mfluxn = rho * velocityV[i,j+1] * (mesh.dx[i])

                    if (mfluxn >= 0)
                        # Geometric Parameters
                        # gs
                        if (j == 1) || (!phi.onoff[i,j-1])
                            gs = 1.0
                        else
                            gs = (mesh.dy[j]) / (mesh.dy[j] + mesh.dy[j-1])
                        end
                        # gn
                        gn = (mesh.dy[j]) / (mesh.dy[j] + mesh.dy[j+1])

                        # UPWIND
                        if (j == 1) || (!phi.onoff[i,j-1])
                            if phi.bounds[i,j]
                                vu = find_bondValue(i, j, id, phi, mesh, bounds, 's')
                            else
                                vu = phi.eval[i,j]
                            end
                        else
                            vu = phi.eval[i,j-1]
                        end

                        # DOWNWIND
                        vd = phi.eval[i,j+1]

                        # CENTER
                        vc = phi.eval[i,j]

                        m_coef = 1.0
                        n_coef = 1.0
                        c_coef = 0.0

                        m_tvd = m_coef * gs
                        n_tvd = n_coef * gn + c_coef

                        cf_vc = (1.0 + 0.5 * m_tvd - 0.5 * n_tvd) * mfluxn
                        cf_vd = (0.5 * n_tvd) * mfluxn
                        cf_vu = (-0.5 * m_tvd) * mfluxn

                        anc = cf_vc
                        and1 = cf_vd
                        and2 = cf_vu

                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i,j+1]
                        AV[n] = and1

                        if (j == 1) || (!phi.onoff[i,j-1])
                            b[id] += -1.0 * and2 * vu
                        else
                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i,j-1]
                            AV[n] = and2
                        end

                    elseif (mfluxn < 0)
                        # Geometric Parameters
                        # gn
                        gn = (mesh.dy[j]) / (mesh.dy[j] + mesh.dy[j+1])
                        # gnn
                        if (j == mesh.m1-1) || (!phi.onoff[i,j+2])
                            gnn = 1.0
                        else
                            gnn = (mesh.dy[j+1]) / (mesh.dy[j+1] + mesh.dy[j+2])
                        end

                        # UPWIND
                        if (j == mesh.m1-1) || (!phi.onoff[i,j+2])
                            if phi.bounds[i,j+1]
                                vu = find_bondValue(i, j+1, phi.gIndex[i,j+1], phi, mesh, bounds, 'n')
                            else
                                vu = phi.eval[i,j+1]
                            end
                        else
                            vu = phi.eval[i,j+2]
                        end

                        # DOWNWIND
                        vd = phi.eval[i,j]

                        # CENTER
                        vc = phi.eval[i,j+1]

                        m_coef = 1.0
                        n_coef = 1.0
                        c_coef = 0.0

                        m_tvd = m_coef * gnn
                        n_tvd = n_coef * (1.0 - gn) + c_coef

                        cf_vc = (1.0 + 0.5 * m_tvd - 0.5 * n_tvd) * mfluxn
                        cf_vd = (0.5 * n_tvd) * mfluxn
                        cf_vu = (-0.5 * m_tvd) * mfluxn

                        anc = cf_vd
                        and1 = cf_vc
                        and2 = cf_vu

                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i,j+1]
                        AV[n] = and1

                        if (j == mesh.m1-1) || (!phi.onoff[i,j+2])
                            b[id] += -1.0 * and2 * vu
                        else
                            n += 1
                            AI[n] = id
                            AJ[n] = phi.gIndex[i,j+2]
                            AV[n] = and2
                        end

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
                        AV[n] = ac + (aec + awc + anc + asc) + 1.0e-20
                    else
                        error("inout = false")
                        # AV[n] = ac + (ae + aee + aw + aww + an + ann + as + ass) - (aew + awe + ans + asn) #ac + (aec + awc + anc + asc)
                    end
                    b[id] += b0
                else
                    n += 1
                    AI[n] = id
                    AJ[n] = id
                    if inout
                        AV[n] = ac + (aec + awc + anc + asc) + 1.0e-20
                    else
                        error("inout = false")
                        # AV[n] = ac + (ae + aee + aw + aww + an + ann + as + ass) - (aew + awe + ans + asn) #ac + (aec + awc + anc + asc)
                    end
                end

                # println("")
                # println("$id")
                # println("Center")
                # println("ac + aec + awc + anc + asc -- b[id]")
                # println("$ac + $aec + $awc + $anc + $asc -- $(b[id])")
                #
                # println("")
                #
                # println("awd1 + awd2 + awd3")
                # println("$awd1 + $awd2 + $awc")
                #
                # println("aed1 + aed2 + aed3")
                # println("$aed1 + $aed2 + $aec")
                #
                # println("asd1 + asd2 + asd3")
                # println("$asd1 + $asd2 + $asc")
                #
                # println("and1 + and2 + and3")
                # println("$and1 + $and2 + $anc")

            end
        end
    end

    A = sparse(AI[1:n], AJ[1:n], AV[1:n])

    return A, b
end
