"""

"""
function _projection_PC_velocity_divergence_ end

function _projection_PC_velocity_divergence_(
    velocity::CSVelocity1D,
    mesh::UnionCSMesh1D,
    deltat::DeltaTime,
    material::CSMaterial1D;
    velocityU::Array{<:AbstractFloat,1} = velocity.fValues.uFace,
    T::Type{<:AbstractFloat} = Float64,
    mthreads::Bool = false,
    transientScheme::Signed = 1,
    interpolation::Signed = 1,
)
    n_equations = maximum_globalIndex(velocity.p)

    bdiv = zeros(T, n_equations)

    if (transientScheme == 1)
        coef = 1.0/deltat.dt1
    elseif (transientScheme == 2)
        coef = (deltat.dt2 / (deltat.dt1 * (deltat.dt1 + deltat.dt2)))
    elseif (transientScheme == 3)
        coef = (1.0 / deltat.dt1 ) + ( 1.0 / (deltat.dt1 + deltat.dt2))
    else
        error("Transient scheme unimplemented...")
    end

    for i in 1:mesh.l1
        if velocity.p.onoff[i,j]
            id = velocity.p.gIndex[i,j]

            # Aux variables
            rhow = 0.0
            rhoe = 0.0
            rhos = 0.0
            rhon = 0.0
            rhob = 0.0
            rhot = 0.0

            # Component x - velocity u
            if (mesh.l1 != 1)
                if (i == 1)
                    rhoe = density_interpolation(
                        mesh.dx[i], mesh.dx[i+1], material.ρ[i], material.ρ[i+1];
                        interpolation = interpolation
                    )
                    rhow = material.ρ[i]

                elseif (i == mesh.l1)
                    rhoe = material.ρ[i]
                    rhow = density_interpolation(
                        mesh.dx[i], mesh.dx[i-1], material.ρ[i], material.ρ[i-1];
                        interpolation = interpolation
                    )

                elseif (!velocity.u.onoff[i-1])
                    rhoe = density_interpolation(
                        mesh.dx[i], mesh.dx[i+1], material.ρ[i], material.ρ[i+1];
                        interpolation = interpolation
                    )
                    rhow = material.ρ[i]

                elseif (!velocity.u.onoff[i+1])
                    rhoe = material.ρ[i]
                    rhow = density_interpolation(
                        mesh.dx[i], mesh.dx[i-1], material.ρ[i], material.ρ[i-1];
                        interpolation = interpolation
                    )

                else
                    rhoe = density_interpolation(
                        mesh.dx[i], mesh.dx[i+1], material.ρ[i], material.ρ[i+1];
                        interpolation = interpolation
                    )
                    rhow = density_interpolation(
                        mesh.dx[i], mesh.dx[i-1], material.ρ[i], material.ρ[i-1];
                        interpolation = interpolation
                    )

                end

                bdiv[id] += coef * (rhow * velocityU[i] - rhoe * velocityU[i+1])

            end

        end
    end

    return bdiv
end

function _projection_PC_velocity_divergence_(
    velocity::CSVelocity2D,
    mesh::UnionCSMesh2D,
    deltat::DeltaTime,
    material::CSMaterial2D;
    velocityU::Array{<:AbstractFloat,2} = velocity.fValues.uFace,
    velocityV::Array{<:AbstractFloat,2} = velocity.fValues.vFace,
    T::Type{<:AbstractFloat} = Float64,
    mthreads::Bool = false,
    transientScheme::Signed = 1,
    interpolation::Signed = 1,
)
    n_equations = maximum_globalIndex(velocity.p)

    bdiv = zeros(T, n_equations)

    if (transientScheme == 1)
        coef = 1.0/deltat.dt1
    elseif (transientScheme == 2)
        coef = (deltat.dt2 / (deltat.dt1 * (deltat.dt1 + deltat.dt2)))
    elseif (transientScheme == 3)
        coef = (1.0 / deltat.dt1 ) + ( 1.0 / (deltat.dt1 + deltat.dt2))
    else
        error("Transient scheme unimplemented...")
    end

    for i in 1:mesh.l1
        for j in 1:mesh.m1
            if velocity.p.onoff[i,j]
                id = velocity.p.gIndex[i,j]

                # Aux variables
                rhow = 0.0
                rhoe = 0.0
                rhos = 0.0
                rhon = 0.0
                rhob = 0.0
                rhot = 0.0

                # Component x - velocity u
                if (mesh.l1 != 1)
                    if (i == 1)
                        rhoe = density_interpolation(
                            mesh.dx[i], mesh.dx[i+1], material.ρ[i,j], material.ρ[i+1,j];
                            interpolation = interpolation
                        )
                        rhow = material.ρ[i,j]

                    elseif (i == mesh.l1)
                        rhoe = material.ρ[i,j]
                        rhow = density_interpolation(
                            mesh.dx[i], mesh.dx[i-1], material.ρ[i,j], material.ρ[i-1,j];
                            interpolation = interpolation
                        )

                    elseif (!velocity.u.onoff[i-1,j])
                        rhoe = density_interpolation(
                            mesh.dx[i], mesh.dx[i+1], material.ρ[i,j], material.ρ[i+1,j];
                            interpolation = interpolation
                        )
                        rhow = material.ρ[i,j]

                    elseif (!velocity.u.onoff[i+1,j])
                        rhoe = material.ρ[i,j]
                        rhow = density_interpolation(
                            mesh.dx[i], mesh.dx[i-1], material.ρ[i,j], material.ρ[i-1,j];
                            interpolation = interpolation
                        )

                    else
                        rhoe = density_interpolation(
                            mesh.dx[i], mesh.dx[i+1], material.ρ[i,j], material.ρ[i+1,j];
                            interpolation = interpolation
                        )
                        rhow = density_interpolation(
                            mesh.dx[i], mesh.dx[i-1], material.ρ[i,j], material.ρ[i-1,j];
                            interpolation = interpolation
                        )

                    end

                    bdiv[id] += coef * (rhow * velocityU[i,j] - rhoe * velocityU[i+1,j]) * mesh.dy[j]

                end

                # Component y - velocity v
                if (mesh.m1 != 1)
                    if (j == 1)
                        rhon = density_interpolation(
                            mesh.dy[j], mesh.dy[j+1], material.ρ[i,j], material.ρ[i,j+1];
                            interpolation = interpolation
                        )
                        rhos = material.ρ[i,j]

                    elseif (j == mesh.m1)
                        rhon = material.ρ[i,j]
                        rhos = density_interpolation(
                            mesh.dy[j], mesh.dy[j-1], material.ρ[i,j], material.ρ[i,j-1];
                            interpolation = interpolation
                        )

                    elseif (!velocity.v.onoff[i,j-1])
                        rhon = density_interpolation(
                            mesh.dy[j], mesh.dy[j+1], material.ρ[i,j], material.ρ[i,j+1];
                            interpolation = interpolation
                        )
                        rhos = material.ρ[i,j]

                    elseif (!velocity.v.onoff[i,j+1])
                        rhon = material.ρ[i,j]
                        rhos = density_interpolation(
                            mesh.dy[j], mesh.dy[j-1], material.ρ[i,j], material.ρ[i,j-1];
                            interpolation = interpolation
                        )

                    else
                        rhon = density_interpolation(
                            mesh.dy[j], mesh.dy[j+1], material.ρ[i,j], material.ρ[i,j+1];
                            interpolation = interpolation
                        )
                        rhos = density_interpolation(
                            mesh.dy[j], mesh.dy[j-1], material.ρ[i,j], material.ρ[i,j-1];
                            interpolation = interpolation
                        )

                    end

                    bdiv[id] += coef * (rhos * velocityV[i,j] - rhon * velocityV[i,j+1]) * mesh.dx[i]

                end

            end
        end
    end

    return bdiv
end

function _projection_PC_velocity_divergence_(
    velocity::CSVelocity3D,
    mesh::UnionCSMesh3D,
    deltat::DeltaTime,
    material::CSMaterial3D;
    velocityU::Array{<:AbstractFloat,3} = velocity.fValues.uFace,
    velocityV::Array{<:AbstractFloat,3} = velocity.fValues.vFace,
    velocityW::Array{<:AbstractFloat,3} = velocity.fValues.wFace,
    T::Type{<:AbstractFloat} = Float64,
    mthreads::Bool = false,
    transientScheme::Signed = 1,
    interpolation::Signed = 1,
)
    n_equations = maximum_globalIndex(velocity.p)

    bdiv = zeros(T, n_equations)

    if (transientScheme == 1)
        coef = 1.0/deltat.dt1
    elseif (transientScheme == 2)
        coef = (deltat.dt2 / (deltat.dt1 * (deltat.dt1 + deltat.dt2)))
    elseif (transientScheme == 3)
        coef = (1.0 / deltat.dt1 ) + ( 1.0 / (deltat.dt1 + deltat.dt2))
    else
        error("Transient scheme unimplemented...")
    end

    for i in 1:mesh.l1
        for j in 1:mesh.m1
            for k in 1:mesh.n1
                if velocity.p.onoff[i,j,k]
                    id = velocity.p.gIndex[i,j,k]

                    # Aux variables
                    rhow = 0.0
                    rhoe = 0.0
                    rhos = 0.0
                    rhon = 0.0
                    rhob = 0.0
                    rhot = 0.0

                    # Component x - velocity u
                    if (mesh.l1 != 1)
                        if (i == 1)
                            rhoe = density_interpolation(
                                mesh.dx[i], mesh.dx[i+1], material.ρ[i,j,k], material.ρ[i+1,j,k];
                                interpolation = interpolation
                            )
                            rhow = material.ρ[i,j,k]

                        elseif (i == mesh.l1)
                            rhoe = material.ρ[i,j,k]
                            rhow = density_interpolation(
                                mesh.dx[i], mesh.dx[i-1], material.ρ[i,j,k], material.ρ[i-1,j,k];
                                interpolation = interpolation
                            )

                        elseif (!velocity.u.onoff[i-1,j,k])
                            rhoe = density_interpolation(
                                mesh.dx[i], mesh.dx[i+1], material.ρ[i,j,k], material.ρ[i+1,j,k];
                                interpolation = interpolation
                            )
                            rhow = material.ρ[i,j,k]

                        elseif (!velocity.u.onoff[i+1,j,k])
                            rhoe = material.ρ[i,j,k]
                            rhow = density_interpolation(
                                mesh.dx[i], mesh.dx[i-1], material.ρ[i,j,k], material.ρ[i-1,j,k];
                                interpolation = interpolation
                            )

                        else
                            rhoe = density_interpolation(
                                mesh.dx[i], mesh.dx[i+1], material.ρ[i,j,k], material.ρ[i+1,j,k];
                                interpolation = interpolation
                            )
                            rhow = density_interpolation(
                                mesh.dx[i], mesh.dx[i-1], material.ρ[i,j,k], material.ρ[i-1,j,k];
                                interpolation = interpolation
                            )

                        end

                        bdiv[id] += coef * (rhow * velocityU[i,j,k] - rhoe * velocityU[i+1,j,k]) * (mesh.dy[j] * mesh.dz[k])

                    end

                    # Component y - velocity v
                    if (mesh.m1 != 1)
                        if (j == 1)
                            rhon = density_interpolation(
                                mesh.dy[j], mesh.dy[j+1], material.ρ[i,j,k], material.ρ[i,j+1,k];
                                interpolation = interpolation
                            )
                            rhos = material.ρ[i,j,k]

                        elseif (j == mesh.m1)
                            rhon = material.ρ[i,j,k]
                            rhos = density_interpolation(
                                mesh.dy[j], mesh.dy[j-1], material.ρ[i,j,k], material.ρ[i,j-1,k];
                                interpolation = interpolation
                            )

                        elseif (!velocity.v.onoff[i,j-1,k])
                            rhon = density_interpolation(
                                mesh.dy[j], mesh.dy[j+1], material.ρ[i,j,k], material.ρ[i,j+1,k];
                                interpolation = interpolation
                            )
                            rhos = material.ρ[i,j,k]

                        elseif (!velocity.v.onoff[i,j+1,k])
                            rhon = material.ρ[i,j,k]
                            rhos = density_interpolation(
                                mesh.dy[j], mesh.dy[j-1], material.ρ[i,j,k], material.ρ[i,j-1,k];
                                interpolation = interpolation
                            )

                        else
                            rhon = density_interpolation(
                                mesh.dy[j], mesh.dy[j+1], material.ρ[i,j,k], material.ρ[i,j+1,k];
                                interpolation = interpolation
                            )
                            rhos = density_interpolation(
                                mesh.dy[j], mesh.dy[j-1], material.ρ[i,j,k], material.ρ[i,j-1,k];
                                interpolation = interpolation
                            )

                        end

                        bdiv[id] += coef * (rhos * velocityV[i,j,k] - rhon * velocityV[i,j+1,k]) * (mesh.dx[i] * mesh.dz[k])

                    end

                    # Component z - velocity w
                    if (mesh.n1 != 1)
                        if (k == 1)
                            rhot = density_interpolation(
                                mesh.dz[k], mesh.dz[k+1], material.ρ[i,j,k], material.ρ[i,j,k+1];
                                interpolation = interpolation
                            )
                            rhob = material.ρ[i,j,k]

                        elseif (k == mesh.n1)
                            rhot = material.ρ[i,j,k]
                            rhob = density_interpolation(
                                mesh.dz[k], mesh.dz[k-1], material.ρ[i,j,k], material.ρ[i,j,k-1];
                                interpolation = interpolation
                            )

                        elseif (!velocity.w.onoff[i,j,k-1])
                            rhot = density_interpolation(
                                mesh.dz[k], mesh.dz[k+1], material.ρ[i,j,k], material.ρ[i,j,k+1];
                                interpolation = interpolation
                            )
                            rhob = material.ρ[i,j,k]

                        elseif (!velocity.w.onoff[i,j,k+1])
                            rhot = material.ρ[i,j,k]
                            rhob = density_interpolation(
                                mesh.dz[k], mesh.dz[k-1], material.ρ[i,j,k], material.ρ[i,j,k-1];
                                interpolation = interpolation
                            )

                        else
                            rhot = density_interpolation(
                                mesh.dz[k], mesh.dz[k+1], material.ρ[i,j,k], material.ρ[i,j,k+1];
                                interpolation = interpolation
                            )
                            rhob = density_interpolation(
                                mesh.dz[k], mesh.dz[k-1], material.ρ[i,j,k], material.ρ[i,j,k-1];
                                interpolation = interpolation
                            )

                        end

                        bdiv[id] += coef * (rhob * velocityW[i,j,k] - rhot * velocityW[i,j,k+1]) * (mesh.dx[i] * mesh.dy[j])

                    end

                end
            end
        end
    end

    return bdiv
end

function _projection_PC_velocity_divergence_(
    velocity::CSVelocity1D,
    mesh::UnionCSMesh1D,
    deltat::DeltaTime,
    material::UnionCSConstantMaterial;
    velocityU::Array{<:AbstractFloat,1} = velocity.fValues.uFace,
    T::Type{<:AbstractFloat} = Float64,
    mthreads::Bool = false,
    transientScheme::Signed = 1,
    interpolation::Signed = 1,
)
    n_equations = maximum_globalIndex(velocity.p)

    bdiv = zeros(T, n_equations)

    if (transientScheme == 1)
        coef = 1.0/deltat.dt1
    elseif (transientScheme == 2)
        coef = (deltat.dt2 / (deltat.dt1 * (deltat.dt1 + deltat.dt2)))
    elseif (transientScheme == 3)
        coef = (1.0 / deltat.dt1 ) + ( 1.0 / (deltat.dt1 + deltat.dt2))
    else
        error("Transient scheme unimplemented...")
    end

    for i in 1:mesh.l1
        if velocity.p.onoff[i,j]
            id = velocity.p.gIndex[i,j]

            if (mesh.l1 != 1)
                bdiv[id] += coef * material.ρ * (velocityU[i+1] - velocityU[i])
            end

        end
    end

    return bdiv
end

function _projection_PC_velocity_divergence_(
    velocity::CSVelocity2D,
    mesh::UnionCSMesh2D,
    deltat::DeltaTime,
    material::UnionCSConstantMaterial;
    velocityU::Array{<:AbstractFloat,2} = velocity.fValues.uFace,
    velocityV::Array{<:AbstractFloat,2} = velocity.fValues.vFace,
    T::Type{<:AbstractFloat} = Float64,
    mthreads::Bool = false,
    transientScheme::Signed = 1,
    interpolation::Signed = 1,
)
    n_equations = maximum_globalIndex(velocity.p)

    bdiv = zeros(T, n_equations)

    if (transientScheme == 1)
        coef = 1.0/deltat.dt1
    elseif (transientScheme == 2)
        coef = (deltat.dt2 / (deltat.dt1 * (deltat.dt1 + deltat.dt2)))
    elseif (transientScheme == 3)
        coef = (1.0 / deltat.dt1 ) + ( 1.0 / (deltat.dt1 + deltat.dt2))
    else
        error("Transient scheme unimplemented...")
    end

    for i in 1:mesh.l1
        for j in 1:mesh.m1
            if velocity.p.onoff[i,j]
                id = velocity.p.gIndex[i,j]

                # Component x - velocity u
                if (mesh.l1 != 1)
                    bdiv[id] += coef * material.ρ * (velocityU[i+1,j] - velocityU[i,j])
                end

                # Component y - velocity v
                if (mesh.m1 != 1)
                    bdiv[id] += coef * material.ρ * (velocityV[i,j+1] - velocityV[i,j])
                end

            end
        end
    end

    return bdiv
end

function _projection_PC_velocity_divergence_(
    velocity::CSVelocity3D,
    mesh::UnionCSMesh3D,
    deltat::DeltaTime,
    material::UnionCSConstantMaterial;
    velocityU::Array{<:AbstractFloat,3} = velocity.fValues.uFace,
    velocityV::Array{<:AbstractFloat,3} = velocity.fValues.vFace,
    velocityW::Array{<:AbstractFloat,3} = velocity.fValues.wFace,
    T::Type{<:AbstractFloat} = Float64,
    mthreads::Bool = false,
    transientScheme::Signed = 1,
    interpolation::Signed = 1,
)
    n_equations = maximum_globalIndex(velocity.p)

    bdiv = zeros(T, n_equations)

    if (transientScheme == 1)
        coef = 1.0/deltat.dt1
    elseif (transientScheme == 2)
        coef = (deltat.dt2 / (deltat.dt1 * (deltat.dt1 + deltat.dt2)))
    elseif (transientScheme == 3)
        coef = (1.0 / deltat.dt1 ) + ( 1.0 / (deltat.dt1 + deltat.dt2))
    else
        error("Transient scheme unimplemented...")
    end

    for i in 1:mesh.l1
        for j in 1:mesh.m1
            for k in 1:mesh.n1
                if velocity.p.onoff[i,j,k]
                    id = velocity.p.gIndex[i,j,k]

                    # Component x - velocity u
                    if (mesh.l1 != 1)
                        bdiv[id] += coef * material.ρ * (velocityU[i+1,j,k] - velocityU[i,j,k])
                    end

                    # Component y - velocity v
                    if (mesh.m1 != 1)
                        bdiv[id] += coef * material.ρ * (velocityV[i,j+1,k] - velocityV[i,j,k])
                    end

                    # Component z - velocity w
                    if (mesh.n1 != 1)
                        bdiv[id] += coef * material.ρ * (velocityW[i,j,k+1] - velocityW[i,j,k])
                    end

                end
            end
        end
    end

    return bdiv
end
