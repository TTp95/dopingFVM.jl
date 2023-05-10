"""

"""
function _discretize_convection_secondorderupwind_TVD_ end


function _discretize_convection_secondorderupwind_TVD_(
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
    # TVD
    # UP
    # m_tvd = 0.0
    # n_tvd = 0.0

    # CD
    # m_tvd = 0.0
    # n_tvd = 1.0

    # SOU
    # m_tvd = 1.0
    # n_tvd = 0.0

    function tvd_coef(rf)
        value_m = 0.0
        value_n = 0.0

        if rf < 0.0
            value_m = 0.0
            value_n = 0.0

        elseif rf > 1.0
            value_m = 0.0
            value_n = 1.0

        else
            value_m = 1.0
            value_n = 0.0

        end

        # UP
        # value_m = 0.0
        # value_n = 0.0

        # DOWNWIND
        # value_m = 0.0
        # value_n = 2.0

        # CD
        # value_m = 0.0
        # value_n = 1.0

        # SOU
        # value_m = 1.0
        # value_n = 0.0

        # FROMM
        # value_m = 0.5
        # value_n = 0.5

        # QUICK
        # value_m = 0.25
        # value_n = 0.75


        return value_m, value_n
    end

    function rf(vu, vc, vd)
        value = 1.0 * (vc - vu)/((vd - vc) + 1.0e-20)

        # println("$(value)")

        return value
    end


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
                b1 = 0.0
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
                mflux = 0.0
                mfluxw = 0.0
                mfluxe = 0.0
                mfluxs = 0.0
                mfluxn = 0.0
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

                    mflux = -1.0 * rho * velocityU[i,j] * (mesh.dy[j])
                    mfluxw = -1.0 * rho * velocityU[i,j] * (mesh.dy[j])

                    if (mflux <= 0)
                        # Velocity values

                        # UPWIND
                        if (i == 2)
                            vu = 0.0 # phi.eval[i-1,j]
                        else
                            vu = phi.eval[i-2,j]
                        end

                        # DOWNWIND
                        vd = phi.eval[i,j]

                        # CENTER
                        vc = phi.eval[i-1,j]

                        vrf = rf(vu, vc, vd)
                        m_tvd, n_tvd = tvd_coef(vrf)

                        # Compute L and K
                        l_tvd = 1.0 + 0.5 * m_tvd - 0.5 * n_tvd
                        k_tvd = 0.5 * n_tvd

                        cf_vc = l_tvd * mflux
                        cf_vd = k_tvd * mflux
                        cf_vu = -0.5 * m_tvd * mflux

                        awc = cf_vd
                        awd1 = cf_vc
                        awd2 = cf_vu

                        b1 = 0.0

                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i-1,j]
                        AV[n] = awd1

                        n += 1
                        AI[n] = id
                        if (i == 2)
                            b[id] += -1.0 * awd2 * vu
                            AJ[n] = phi.gIndex[i-1,j]
                            AV[n] = 0.0 # awd2
                        else
                            AJ[n] = phi.gIndex[i-2,j]
                            AV[n] = awd2
                        end

                        b[id] += b1

                    end

                    if (mflux > 0)
                        # Velocity values

                        # UPWIND
                        if (i == mesh.l1)
                            vu = 0.0 # phi.eval[i,j]
                        else
                            vu = phi.eval[i+1,j]
                        end

                        # DOWNWIND
                        vd = phi.eval[i-1,j]

                        # CENTER
                        vc = phi.eval[i,j]

                        vrf = rf(vu, vc, vd)
                        m_tvd, n_tvd = tvd_coef(vrf)

                        # Compute L and K
                        l_tvd = 1.0 + 0.5 * m_tvd - 0.5 * n_tvd
                        k_tvd = 0.5 * n_tvd

                        cf_vc = l_tvd * mflux
                        cf_vd = k_tvd * mflux
                        cf_vu = -0.5 * m_tvd * mflux

                        awc = cf_vc
                        awd1 = cf_vd
                        awd2 = cf_vu

                        b1 = 0.0

                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i-1,j]
                        AV[n] = awd1

                        n += 1
                        AI[n] = id
                        if (i == mesh.l1)
                            b[id] += -1.0 * awd2 * vu
                            AJ[n] = phi.gIndex[i,j]
                            AV[n] = 0.0 # awd2
                        else
                            AJ[n] = phi.gIndex[i+1,j]
                            AV[n] = awd2
                        end

                        b[id] += b1

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

                    if (mflux >= 0)
                        # Velocity values

                        # UPWIND
                        if (i == 1)
                            vu = 0.0 # phi.eval[i,j]
                        else
                            vu = phi.eval[i-1,j]
                        end

                        # DOWNWIND
                        vd = phi.eval[i+1,j]

                        # CENTER
                        vc = phi.eval[i,j]

                        vrf = rf(vu, vc, vd)
                        m_tvd, n_tvd = tvd_coef(vrf)

                        # Compute L and K
                        l_tvd = 1.0 + 0.5 * m_tvd - 0.5 * n_tvd
                        k_tvd = 0.5 * n_tvd

                        cf_vc = l_tvd * mflux
                        cf_vd = k_tvd * mflux
                        cf_vu = -0.5 * m_tvd * mflux

                        aec = cf_vc
                        aed1 = cf_vd
                        aed2 = cf_vu

                        b1 = 0.0

                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i+1,j]
                        AV[n] = aed1

                        n += 1
                        AI[n] = id
                        if (i == 1)
                            b[id] += -1.0 * aed2 * vu
                            AJ[n] = phi.gIndex[i,j]
                            AV[n] = 0.0 # aed2
                        else
                            AJ[n] = phi.gIndex[i-1,j]
                            AV[n] = aed2
                        end

                        b[id] += b1

                    end

                    if (mflux < 0)
                        # Velocity values

                        # UPWIND
                        if (i == mesh.l1-1)
                            vu = 0.0 # phi.eval[i+1,j]
                        else
                            vu = phi.eval[i+2,j]
                        end

                        # DOWNWIND
                        vd = phi.eval[i,j]

                        # CENTER
                        vc = phi.eval[i+1,j]

                        vrf = rf(vu, vc, vd)
                        m_tvd, n_tvd = tvd_coef(vrf)

                        # Compute L and K
                        l_tvd = 1.0 + 0.5 * m_tvd - 0.5 * n_tvd
                        k_tvd = 0.5 * n_tvd

                        cf_vc = l_tvd * mflux
                        cf_vd = k_tvd * mflux
                        cf_vu = -0.5 * m_tvd * mflux

                        aec = cf_vd
                        aed1 = cf_vc
                        aed2 = cf_vu

                        b1 = 0.0

                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i+1,j]
                        AV[n] = aed1

                        n += 1
                        AI[n] = id
                        if (i == mesh.l1-1)
                            b[id] += -1.0 * aed2 * vu
                            AJ[n] = phi.gIndex[i+1,j]
                            AV[n] = 0.0 # aed2
                        else
                            AJ[n] = phi.gIndex[i+2,j]
                            AV[n] = aed2
                        end

                        b[id] += b1

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


                    if (mflux <= 0)
                        # Velocity values

                        # UPWIND
                        if (j == 2)
                            vu = 0.0 # phi.eval[i,j-1]
                        else
                            vu = phi.eval[i,j-2]
                        end

                        # DOWNWIND
                        vd = phi.eval[i,j]

                        # CENTER
                        vc = phi.eval[i,j-1]

                        vrf = rf(vu, vc, vd)
                        m_tvd, n_tvd = tvd_coef(vrf)

                        # Compute L and K
                        l_tvd = 1.0 + 0.5 * m_tvd - 0.5 * n_tvd
                        k_tvd = 0.5 * n_tvd

                        cf_vc = l_tvd * mflux
                        cf_vd = k_tvd * mflux
                        cf_vu = -0.5 * m_tvd * mflux

                        asc = cf_vd
                        asd1 = cf_vc
                        asd2 = cf_vu

                        b1 = 0.0

                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i,j-1]
                        AV[n] = asd1

                        n += 1
                        AI[n] = id
                        if (j == 2)
                            b[id] += -1.0 * asd2 * vu
                            AJ[n] = phi.gIndex[i,j-1]
                            AV[n] = 0.0 # asd2
                        else
                            AJ[n] = phi.gIndex[i,j-2]
                            AV[n] = asd2
                        end

                        b[id] += b1

                    end

                    if (mflux > 0)
                        # Velocity values

                        # UPWIND
                        if (j == mesh.m1)
                            vu = 0.0 # phi.eval[i,j]
                        else
                            vu = phi.eval[i,j+1]
                        end

                        # DOWNWIND
                        vd = phi.eval[i,j-1]

                        # CENTER
                        vc = phi.eval[i,j]

                        vrf = rf(vu, vc, vd)
                        m_tvd, n_tvd = tvd_coef(vrf)

                        # Compute L and K
                        l_tvd = 1.0 + 0.5 * m_tvd - 0.5 * n_tvd
                        k_tvd = 0.5 * n_tvd

                        cf_vc = l_tvd * mflux
                        cf_vd = k_tvd * mflux
                        cf_vu = -0.5 * m_tvd * mflux

                        asc = cf_vc
                        asd1 = cf_vd
                        asd2 = cf_vu

                        b1 = 0.0

                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i,j-1]
                        AV[n] = asd1

                        n += 1
                        AI[n] = id
                        if (j == mesh.m1)
                            b[id] += -1.0 * asd2 * vu
                            AJ[n] = phi.gIndex[i,j]
                            AV[n] = 0.0 # asd2
                        else
                            AJ[n] = phi.gIndex[i,j+1]
                            AV[n] = asd2
                        end

                        b[id] += b1

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

                    if (mflux >= 0)
                        # Velocity values

                        # UPWIND
                        if (j == 1)
                            vu = 0.0 # phi.eval[i,j]
                        else
                            vu = phi.eval[i,j-1]
                        end

                        # DOWNWIND
                        vd = phi.eval[i,j+1]

                        # CENTER
                        vc = phi.eval[i,j]

                        vrf = rf(vu, vc, vd)
                        m_tvd, n_tvd = tvd_coef(vrf)

                        # Compute L and K
                        l_tvd = 1.0 + 0.5 * m_tvd - 0.5 * n_tvd
                        k_tvd = 0.5 * n_tvd

                        cf_vc = l_tvd * mflux
                        cf_vd = k_tvd * mflux
                        cf_vu = -0.5 * m_tvd * mflux

                        anc = cf_vc
                        and1 = cf_vd
                        and2 = cf_vu

                        b1 = 0.0

                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i,j+1]
                        AV[n] = and1

                        n += 1
                        AI[n] = id
                        if (j == 1)
                            b[id] += -1.0 * and2 * vu
                            AJ[n] = phi.gIndex[i,j]
                            AV[n] = 0.0 # and2
                        else
                            AJ[n] = phi.gIndex[i,j-1]
                            AV[n] = and2
                        end

                        b[id] += b1

                    end

                    if (mflux < 0)
                        # Velocity values

                        # UPWIND
                        if (j == mesh.m1-1)
                            vu = 0.0 # phi.eval[i,j+1]
                        else
                            vu = phi.eval[i,j+2]
                        end

                        # DOWNWIND
                        vd = phi.eval[i,j]

                        # CENTER
                        vc = phi.eval[i,j+1]

                        vrf = rf(vu, vc, vd)
                        m_tvd, n_tvd = tvd_coef(vrf)

                        # Compute L and K
                        l_tvd = 1.0 + 0.5 * m_tvd - 0.5 * n_tvd
                        k_tvd = 0.5 * n_tvd

                        cf_vc = l_tvd * mflux
                        cf_vd = k_tvd * mflux
                        cf_vu = -0.5 * m_tvd * mflux

                        anc = cf_vd
                        and1 = cf_vc
                        and2 = cf_vu

                        b1 = 0.0

                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i,j+1]
                        AV[n] = and1

                        n += 1
                        AI[n] = id
                        if (j == mesh.m1-1)
                            b[id] += -1.0 * and2 * vu
                            AJ[n] = phi.gIndex[i,j+1]
                            AV[n] = 0.0 # and2
                        else
                            AJ[n] = phi.gIndex[i,j+2]
                            AV[n] = and2
                        end

                        b[id] += b1

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
                        # AV[n] = ac + (aec + awc + anc + asc) + 1.0e-20


                        #AV[n] = ac + (ae + aee + aw + aww + an + ann + as + ass) - (aew + awe + ans + asn) + (acfw + acfe + acfs + acfn) + (mfluxw + mfluxe + mfluxs + mfluxn)
                        # println("$(AV[n]) -> $i, $j -> $ac + $aec + $awc + $anc + $asc")
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
                        # AV[n] = ac + (aec + awc + anc + asc) + 1.0e-20


                        #AV[n] = ac + (ae + aee + aw + aww + an + ann + as + ass) - (aew + awe + ans + asn) + (acfw + acfe + acfs + acfn) + (mfluxw + mfluxe + mfluxs + mfluxn)
                        # println("$(AV[n]) -> $i, $j -> $ac + $aec + $awc + $anc + $asc")
                    else
                        error("inout = false")
                        # AV[n] = ac + (ae + aee + aw + aww + an + ann + as + ass) - (aew + awe + ans + asn) #ac + (aec + awc + anc + asc)
                    end
                end
            end
        end
    end

    A = sparse(AI[1:n], AJ[1:n], AV[1:n])

    return A, b
end
