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
    # TVD functions
    lim(r) = r
    rf(vu, vc, vd) = (vc - vu)/(vd - vc)

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
                # aw = 0.0
                # awe = 0.0
                # aww = 0.0
                awc = 0.0
                awd = 0.0
                # b1 = 0.0
                # ae = 0.0
                # aew = 0.0
                # aee = 0.0
                aec = 0.0
                aed = 0.0
                # b2 = 0.0
                # as = 0.0
                # asn = 0.0
                # ass = 0.0
                asc = 0.0
                asd = 0.0
                # b3 = 0.0
                # an = 0.0
                # ans = 0.0
                # ann = 0.0
                anc = 0.0
                and = 0.0
                # b4 = 0.0
                mflux = 0.0
                mfluxw = 0.0
                mfluxe = 0.0
                mfluxs = 0.0
                mfluxn = 0.0
                # acfw = 0.0
                # acfe = 0.0
                # acfs = 0.0
                # acfn = 0.0
                # num = 0.0
                # den = 0.0
                vu = 0.0
                vd = 0.0
                vc = 0.0

                #West Coefficents
                if (i != 1) && (mesh.l1 != 1) && (phi.onoff[i-1,j])
                    rho = density_interpolation(
                        mesh.dx[i], mesh.dx[i-1], material.ρ[i,j], material.ρ[i-1,j];
                        interpolation = interpolation
                    )

                    mflux = -1.0 * rho * velocityU[i,j] * (mesh.dy[j])
                    mfluxw = -1.0 * rho * velocityU[i,j] * (mesh.dy[j])

                    if (mflux >= 0)
                        # Velocity values

                        # UPWIND
                        if (i == 2)
                            vu = phi.eval[i-1,j]
                        else
                            vu = phi.eval[i-2,j]
                        end

                        # DOWNWIND
                        vd = phi.eval[i,j]

                        # CENTER
                        vc = phi.eval[i-1,j]

                        vrf = rf(vu, vd, vc)
                        vlim = lim(vrf)

                        awc = (0.5 * vlim) * mflux
                        awd = (1.0 - 0.5 * vlim) * mflux
                        b1 = 0.0

                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i-1,j]
                        AV[n] = awd
                        b[id] += b1

                    end

                    if (mflux < 0)
                        # Velocity values

                        # UPWIND
                        if (i == mesh.l1)
                            vu = phi.eval[i,j]
                        else
                            vu = phi.eval[i+1,j]
                        end

                        # DOWNWIND
                        vd = phi.eval[i-1,j]

                        # CENTER
                        vc = phi.eval[i,j]

                        vrf = rf(vu, vd, vc)
                        vlim = lim(vrf)

                        awc = (1.0 - 0.5 * vlim) * mflux
                        awd = (0.5 * vlim) * mflux
                        b1 = 0.0

                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i-1,j]
                        AV[n] = awd
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
                            vu = phi.eval[i,j]
                        else
                            vu = phi.eval[i-1,j]
                        end

                        # DOWNWIND
                        vd = phi.eval[i+1,j]

                        # CENTER
                        vc = phi.eval[i,j]

                        vrf = rf(vu, vd, vc)
                        vlim = lim(vrf)

                        aec = (1.0 - 0.5 * vlim) * mflux
                        aed = (0.5 * vlim) * mflux
                        b1 = 0.0

                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i+1,j]
                        AV[n] = aed
                        b[id] += b1

                    end

                    if (mflux < 0)
                        # Velocity values

                        # UPWIND
                        if (i == mesh.l1-1)
                            vu = phi.eval[i+1,j]
                        else
                            vu = phi.eval[i+2,j]
                        end

                        # DOWNWIND
                        vd = phi.eval[i,j]

                        # CENTER
                        vc = phi.eval[i+1,j]

                        vrf = rf(vu, vd, vc)
                        vlim = lim(vrf)

                        aec = (0.5 * vlim) * mflux
                        aed = (1.0 - 0.5 * vlim) * mflux
                        b1 = 0.0

                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i+1,j]
                        AV[n] = aed
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

                    if (mflux >= 0)
                        # Velocity values

                        # UPWIND
                        if (j == 2)
                            vu = phi.eval[i,j-1]
                        else
                            vu = phi.eval[i,j-2]
                        end

                        # DOWNWIND
                        vd = phi.eval[i,j]

                        # CENTER
                        vc = phi.eval[i,j-1]

                        vrf = rf(vu, vd, vc)
                        vlim = lim(vrf)

                        asc = (0.5 * vlim) * mflux
                        asd = (1.0 - 0.5 * vlim) * mflux
                        b1 = 0.0

                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i,j-1]
                        AV[n] = asd
                        b[id] += b1

                    end

                    if (mflux < 0)
                        # Velocity values

                        # UPWIND
                        if (j == mesh.m1)
                            vu = phi.eval[i,j]
                        else
                            vu = phi.eval[i,j+1]
                        end

                        # DOWNWIND
                        vd = phi.eval[i,j-1]

                        # CENTER
                        vc = phi.eval[i,j]

                        vrf = rf(vu, vd, vc)
                        vlim = lim(vrf)

                        asc = (1.0 - 0.5 * vlim) * mflux
                        asd = (0.5 * vlim) * mflux
                        b1 = 0.0

                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i,j-1]
                        AV[n] = asd
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
                            vu = phi.eval[i,j]
                        else
                            vu = phi.eval[i,j-1]
                        end

                        # DOWNWIND
                        vd = phi.eval[i,j+1]

                        # CENTER
                        vc = phi.eval[i,j]

                        vrf = rf(vu, vd, vc)
                        vlim = lim(vrf)

                        anc = (1.0 - 0.5 * vlim) * mflux
                        and = (0.5 * vlim) * mflux
                        b1 = 0.0

                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i,j+1]
                        AV[n] = and
                        b[id] += b1

                    end

                    if (mflux < 0)
                        # Velocity values

                        # UPWIND
                        if (j == mesh.m1-1)
                            vu = phi.eval[i,j+1]
                        else
                            vu = phi.eval[i,j+2]
                        end

                        # DOWNWIND
                        vd = phi.eval[i,j]

                        # CENTER
                        vc = phi.eval[i,j+1]

                        vrf = rf(vu, vd, vc)
                        vlim = lim(vrf)

                        anc = (0.5 * vlim) * mflux
                        and = (1.0 - 0.5 * vlim) * mflux
                        b1 = 0.0

                        n += 1
                        AI[n] = id
                        AJ[n] = phi.gIndex[i,j+1]
                        AV[n] = and
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
                        AV[n] = ac + (aec + awc + anc + asc)
                        #AV[n] = ac + (ae + aee + aw + aww + an + ann + as + ass) - (aew + awe + ans + asn) + (acfw + acfe + acfs + acfn) + (mfluxw + mfluxe + mfluxs + mfluxn)
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
                        AV[n] = ac + (aec + awc + anc + asc)
                        #AV[n] = ac + (ae + aee + aw + aww + an + ann + as + ass) - (aew + awe + ans + asn) + (acfw + acfe + acfs + acfn) + (mfluxw + mfluxe + mfluxs + mfluxn)
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
