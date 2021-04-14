"""

"""
function divergence_velocityToArray end

function divergence_velocityToArray(
    velocity::CSVelocity1D,
    mesh::UnionCSMesh1D,
    material::CSMaterial1D;
    velocityU::Array{<:AbstractFloat,1} = velocity.fValues.uFace,
    T::Type{<:AbstractFloat} = Float64,
    mthreads::Bool = false,
    interpolation::Signed = 2,
)
    divArray = zeros(T, mesh.l1, mesh.m1, mesh.n1)

    for i in 1:mesh.l1
        if velocity.p.onoff[i,j]
            id = velocity.p.gIndex[i,j]

            # Aux variables
            gammaw = 0.0
            gammae = 0.0
            gammas = 0.0
            gamman = 0.0
            gammab = 0.0
            gammat = 0.0

            # Component x - velocity u
            if (mesh.l1 != 1)
                if (i == 1)
                    gammae = gamma_interpolation(
                        mesh.dx[i], mesh.dx[i+1], material.Γ[i], material.Γ[i+1];
                        interpolation = interpolation
                    )
                    gammaw = material.Γ[i]

                elseif (i == mesh.l1)
                    gammae = material.Γ[i]
                    gammaw = gamma_interpolation(
                        mesh.dx[i], mesh.dx[i-1], material.Γ[i], material.Γ[i-1];
                        interpolation = interpolation
                    )

                elseif (!velocity.u.onoff[i-1])
                    gammae = gamma_interpolation(
                        mesh.dx[i], mesh.dx[i+1], material.Γ[i], material.Γ[i+1];
                        interpolation = interpolation
                    )
                    gammaw = material.Γ[i]

                elseif (!velocity.u.onoff[i+1])
                    gammae = material.Γ[i]
                    gammaw = gamma_interpolation(
                        mesh.dx[i], mesh.dx[i-1], material.Γ[i], material.Γ[i-1];
                        interpolation = interpolation
                    )

                else
                    gammae = gamma_interpolation(
                        mesh.dx[i], mesh.dx[i+1], material.Γ[i], material.Γ[i+1];
                        interpolation = interpolation
                    )
                    gammaw = gamma_interpolation(
                        mesh.dx[i], mesh.dx[i-1], material.Γ[i], material.Γ[i-1];
                        interpolation = interpolation
                    )

                end

                divArray[i,j] += (gammae * velocityU[i+1] - gammaw * velocityU[i]) / mesh.dx[i]

            end

        end
    end

    return divArray
end

function divergence_velocityToArray(
    velocity::CSVelocity2D,
    mesh::UnionCSMesh2D,
    material::CSMaterial2D;
    velocityU::Array{<:AbstractFloat,2} = velocity.fValues.uFace,
    velocityV::Array{<:AbstractFloat,2} = velocity.fValues.vFace,
    T::Type{<:AbstractFloat} = Float64,
    mthreads::Bool = false,
    interpolation::Signed = 2,
)
    divArray = zeros(T, mesh.l1, mesh.m1)

    for i in 1:mesh.l1
        for j in 1:mesh.m1
            if velocity.p.onoff[i,j]
                id = velocity.p.gIndex[i,j]

                # Aux variables
                gammaw = 0.0
                gammae = 0.0
                gammas = 0.0
                gamman = 0.0
                gammab = 0.0
                gammat = 0.0

                # Component x - velocity u
                if (mesh.l1 != 1)
                    if (i == 1)
                        gammae = gamma_interpolation(
                            mesh.dx[i], mesh.dx[i+1], material.Γ[i,j], material.Γ[i+1,j];
                            interpolation = interpolation
                        )
                        gammaw = material.Γ[i,j]

                    elseif (i == mesh.l1)
                        gammae = material.Γ[i,j]
                        gammaw = gamma_interpolation(
                            mesh.dx[i], mesh.dx[i-1], material.Γ[i,j], material.Γ[i-1,j];
                            interpolation = interpolation
                        )

                    elseif (!velocity.u.onoff[i-1,j])
                        gammae = gamma_interpolation(
                            mesh.dx[i], mesh.dx[i+1], material.Γ[i,j], material.Γ[i+1,j];
                            interpolation = interpolation
                        )
                        gammaw = material.Γ[i,j]

                    elseif (!velocity.u.onoff[i+1,j])
                        gammae = material.Γ[i,j]
                        gammaw = gamma_interpolation(
                            mesh.dx[i], mesh.dx[i-1], material.Γ[i,j], material.Γ[i-1,j];
                            interpolation = interpolation
                        )

                    else
                        gammae = gamma_interpolation(
                            mesh.dx[i], mesh.dx[i+1], material.Γ[i,j], material.Γ[i+1,j];
                            interpolation = interpolation
                        )
                        gammaw = gamma_interpolation(
                            mesh.dx[i], mesh.dx[i-1], material.Γ[i,j], material.Γ[i-1,j];
                            interpolation = interpolation
                        )

                    end

                    divArray[i,j] += (gammae * velocityU[i+1,j] - gammaw * velocityU[i,j]) / mesh.dx[i]

                end

                # Component y - velocity v
                if (mesh.m1 != 1)
                    if (j == 1)
                        gamman = gamma_interpolation(
                            mesh.dy[j], mesh.dy[j+1], material.Γ[i,j], material.Γ[i,j+1];
                            interpolation = interpolation
                        )
                        gammas = material.Γ[i,j]

                    elseif (j == mesh.m1)
                        gamman = material.Γ[i,j]
                        gammas = gamma_interpolation(
                            mesh.dy[j], mesh.dy[j-1], material.Γ[i,j], material.Γ[i,j-1];
                            interpolation = interpolation
                        )

                    elseif (!velocity.v.onoff[i,j-1])
                        gamman = gamma_interpolation(
                            mesh.dy[j], mesh.dy[j+1], material.Γ[i,j], material.Γ[i,j+1];
                            interpolation = interpolation
                        )
                        gammas = material.Γ[i,j]

                    elseif (!velocity.v.onoff[i,j+1])
                        gamman = material.Γ[i,j]
                        gammas = gamma_interpolation(
                            mesh.dy[j], mesh.dy[j-1], material.Γ[i,j], material.Γ[i,j-1];
                            interpolation = interpolation
                        )

                    else
                        gamman = gamma_interpolation(
                            mesh.dy[j], mesh.dy[j+1], material.Γ[i,j], material.Γ[i,j+1];
                            interpolation = interpolation
                        )
                        gammas = gamma_interpolation(
                            mesh.dy[j], mesh.dy[j-1], material.Γ[i,j], material.Γ[i,j-1];
                            interpolation = interpolation
                        )

                    end

                    divArray[i,j] += (gamman * velocityV[i,j+1] - gammas * velocityV[i,j]) / mesh.dy[j]

                end

            end
        end
    end

    return divArray
end

function divergence_velocityToArray(
    velocity::CSVelocity3D,
    mesh::UnionCSMesh3D,
    material::CSMaterial3D;
    velocityU::Array{<:AbstractFloat,3} = velocity.fValues.uFace,
    velocityV::Array{<:AbstractFloat,3} = velocity.fValues.vFace,
    velocityW::Array{<:AbstractFloat,3} = velocity.fValues.wFace,
    T::Type{<:AbstractFloat} = Float64,
    mthreads::Bool = false,
    interpolation::Signed = 2,
)
    divArray = zeros(T, mesh.l1, mesh.m1, mesh.n1)

    for i in 1:mesh.l1
        for j in 1:mesh.m1
            for k in 1:mesh.n1
                if velocity.p.onoff[i,j,k]
                    id = velocity.p.gIndex[i,j,k]

                    # Aux variables
                    gammaw = 0.0
                    gammae = 0.0
                    gammas = 0.0
                    gamman = 0.0
                    gammab = 0.0
                    gammat = 0.0

                    # Component x - velocity u
                    if (mesh.l1 != 1)
                        if (i == 1)
                            gammae = gamma_interpolation(
                                mesh.dx[i], mesh.dx[i+1], material.Γ[i,j,k], material.Γ[i+1,j,k];
                                interpolation = interpolation
                            )
                            gammaw = material.Γ[i,j,k]

                        elseif (i == mesh.l1)
                            gammae = material.Γ[i,j,k]
                            gammaw = gamma_interpolation(
                                mesh.dx[i], mesh.dx[i-1], material.Γ[i,j,k], material.Γ[i-1,j,k];
                                interpolation = interpolation
                            )

                        elseif (!velocity.u.onoff[i-1,j,k])
                            gammae = gamma_interpolation(
                                mesh.dx[i], mesh.dx[i+1], material.Γ[i,j,k], material.Γ[i+1,j,k];
                                interpolation = interpolation
                            )
                            gammaw = material.Γ[i,j,k]

                        elseif (!velocity.u.onoff[i+1,j,k])
                            gammae = material.Γ[i,j,k]
                            gammaw = gamma_interpolation(
                                mesh.dx[i], mesh.dx[i-1], material.Γ[i,j,k], material.Γ[i-1,j,k];
                                interpolation = interpolation
                            )

                        else
                            gammae = gamma_interpolation(
                                mesh.dx[i], mesh.dx[i+1], material.Γ[i,j,k], material.Γ[i+1,j,k];
                                interpolation = interpolation
                            )
                            gammaw = gamma_interpolation(
                                mesh.dx[i], mesh.dx[i-1], material.Γ[i,j,k], material.Γ[i-1,j,k];
                                interpolation = interpolation
                            )

                        end

                        divArray[i,j] += (gammae * velocityU[i+1,j,k] - gammaw * velocityU[i,j,k]) / mesh.dx[i]

                    end

                    # Component y - velocity v
                    if (mesh.m1 != 1)
                        if (j == 1)
                            gamman = gamma_interpolation(
                                mesh.dy[j], mesh.dy[j+1], material.Γ[i,j,k], material.Γ[i,j+1,k];
                                interpolation = interpolation
                            )
                            gammas = material.Γ[i,j,k]

                        elseif (j == mesh.m1)
                            gamman = material.Γ[i,j,k]
                            gammas = gamma_interpolation(
                                mesh.dy[j], mesh.dy[j-1], material.Γ[i,j,k], material.Γ[i,j-1,k];
                                interpolation = interpolation
                            )

                        elseif (!velocity.v.onoff[i,j-1,k])
                            gamman = gamma_interpolation(
                                mesh.dy[j], mesh.dy[j+1], material.Γ[i,j,k], material.Γ[i,j+1,k];
                                interpolation = interpolation
                            )
                            gammas = material.Γ[i,j,k]

                        elseif (!velocity.v.onoff[i,j+1,k])
                            gamman = material.Γ[i,j,k]
                            gammas = gamma_interpolation(
                                mesh.dy[j], mesh.dy[j-1], material.Γ[i,j,k], material.Γ[i,j-1,k];
                                interpolation = interpolation
                            )

                        else
                            gamman = gamma_interpolation(
                                mesh.dy[j], mesh.dy[j+1], material.Γ[i,j,k], material.Γ[i,j+1,k];
                                interpolation = interpolation
                            )
                            gammas = gamma_interpolation(
                                mesh.dy[j], mesh.dy[j-1], material.Γ[i,j,k], material.Γ[i,j-1,k];
                                interpolation = interpolation
                            )

                        end

                        divArray[i,j] += (gamman * velocityV[i,j+1,k] - gammas * velocityV[i,j,k]) / mesh.dy[j]

                    end

                    # Component z - velocity w
                    if (mesh.n1 != 1)
                        if (k == 1)
                            gammat = gamma_interpolation(
                                mesh.dz[k], mesh.dz[k+1], material.Γ[i,j,k], material.Γ[i,j,k+1];
                                interpolation = interpolation
                            )
                            gammab = material.Γ[i,j,k]

                        elseif (k == mesh.n1)
                            gammat = material.Γ[i,j,k]
                            gammab = gamma_interpolation(
                                mesh.dz[k], mesh.dz[k-1], material.Γ[i,j,k], material.Γ[i,j,k-1];
                                interpolation = interpolation
                            )

                        elseif (!velocity.w.onoff[i,j,k-1])
                            gammat = gamma_interpolation(
                                mesh.dz[k], mesh.dz[k+1], material.Γ[i,j,k], material.Γ[i,j,k+1];
                                interpolation = interpolation
                            )
                            gammab = material.Γ[i,j,k]

                        elseif (!velocity.w.onoff[i,j,k+1])
                            gammat = material.Γ[i,j,k]
                            gammab = gamma_interpolation(
                                mesh.dz[k], mesh.dz[k-1], material.Γ[i,j,k], material.Γ[i,j,k-1];
                                interpolation = interpolation
                            )

                        else
                            gammat = gamma_interpolation(
                                mesh.dz[k], mesh.dz[k+1], material.Γ[i,j,k], material.Γ[i,j,k+1];
                                interpolation = interpolation
                            )
                            gammab = gamma_interpolation(
                                mesh.dz[k], mesh.dz[k-1], material.Γ[i,j,k], material.Γ[i,j,k-1];
                                interpolation = interpolation
                            )

                        end

                        divArray[i,j] += (gammat * velocityW[i,j,k+1] - gammab * velocityW[i,j,k]) / mesh.dz[k]

                    end

                end
            end
        end
    end

    return divArray
end

function divergence_velocityToArray(
    velocity::CSVelocity1D,
    mesh::UnionCSMesh1D,
    material::UnionCSConstantMaterial;
    velocityU::Array{<:AbstractFloat,1} = velocity.fValues.uFace,
    T::Type{<:AbstractFloat} = Float64,
    mthreads::Bool = false,
    interpolation::Signed = 2,
)
    divArray = zeros(T, mesh.l1)

    for i in 1:mesh.l1
        if velocity.p.onoff[i,j]
            id = velocity.p.gIndex[i,j]

            if (mesh.l1 != 1)
                divArray[i,j] += material.Γ * (velocityU[i+1] - velocityU[i]) / mesh.dx[i]
            end

        end
    end

    return divArray
end

function divergence_velocityToArray(
    velocity::CSVelocity2D,
    mesh::UnionCSMesh2D,
    material::UnionCSConstantMaterial;
    velocityU::Array{<:AbstractFloat,2} = velocity.fValues.uFace,
    velocityV::Array{<:AbstractFloat,2} = velocity.fValues.vFace,
    T::Type{<:AbstractFloat} = Float64,
    mthreads::Bool = false,
    interpolation::Signed = 2,
)
    divArray = zeros(T, mesh.l1, mesh.m1)

    for i in 1:mesh.l1
        for j in 1:mesh.m1
            if velocity.p.onoff[i,j]
                id = velocity.p.gIndex[i,j]

                # Component x - velocity u
                if (mesh.l1 != 1)
                    divArray[i,j] += material.Γ * (velocityU[i+1,j] - velocityU[i,j]) / mesh.dx[i]
                end

                # Component y - velocity v
                if (mesh.m1 != 1)
                    divArray[i,j] += material.Γ * (velocityV[i,j+1] - velocityV[i,j]) / mesh.dy[j]
                end

            end
        end
    end

    return divArray
end

function divergence_velocityToArray(
    velocity::CSVelocity3D,
    mesh::UnionCSMesh3D,
    material::UnionCSConstantMaterial;
    velocityU::Array{<:AbstractFloat,3} = velocity.fValues.uFace,
    velocityV::Array{<:AbstractFloat,3} = velocity.fValues.vFace,
    velocityW::Array{<:AbstractFloat,3} = velocity.fValues.wFace,
    T::Type{<:AbstractFloat} = Float64,
    mthreads::Bool = false,
    interpolation::Signed = 2,
)
    divArray = zeros(T, mesh.l1, mesh.m1, mesh.n1)

    for i in 1:mesh.l1
        for j in 1:mesh.m1
            for k in 1:mesh.n1
                if velocity.p.onoff[i,j,k]

                    # Component x - velocity u
                    if (mesh.l1 != 1)
                        divArray[i,j] += material.Γ * (velocityU[i+1,j,k] - velocityU[i,j,k]) / mesh.dx[i]
                    end

                    # Component y - velocity v
                    if (mesh.m1 != 1)
                        divArray[i,j] += material.Γ * (velocityV[i,j+1,k] - velocityV[i,j,k]) / mesh.dy[j]
                    end

                    # Component z - velocity w
                    if (mesh.n1 != 1)
                        divArray[i,j] += material.Γ * (velocityW[i,j,k+1] - velocityW[i,j,k]) / mesh.dz[k]
                    end

                end
            end
        end
    end

    return divArray
end
