"""

"""
function PISO_pseudovelocity end

# function PISO_pseudovelocity(
    # AU::Union{
    # SparseVector{<:AbstractFloat,<:Signed},
    # SparseMatrixCSC{<:AbstractFloat,<:Signed},
    # Array{<:AbstractFloat,2},
    # },
    # bU::Array{<:AbstractFloat,1},
    # velocity::CSVelocity1D,
    # mesh::UnionCSMesh1D;
    # relaxP = 1.0,
    # T::Type{<:AbstractFloat} = Float64,
    # full::Bool = true,
# )
    # pseudo_u = zeros(mesh.l1)

    # u = phi_to_vector(velocity.u, mesh)

    # asum_u = 0.0

    # for i in 1:mesh.l1
        # if (velocity.u.onoff[i] == 1.0)

            # id = velocity.u.gIndex[i]

            # asum_u = 0.0

            # if i != 1
                # if (velocity.u.onoff[i-1] == 1.0)
                    # asum_u += AU[id, velocity.u.gIndex[i-1]] * u[velocity.u.gIndex[i-1]]
                # end
            # end

            # if i != mesh.l1
                # if (velocity.u.onoff[i+1] == 1.0)
                    # asum_u += AU[id, velocity.u.gIndex[i+1]] * u[velocity.u.gIndex[i+1]]
                # end
            # end

            # if full
                # if i > 2
                    # if (velocity.u.onoff[i-2] == 1.0)
                        # asum_u += AU[id, velocity.u.gIndex[i-2]] * u[velocity.u.gIndex[i-2]]
                    # end
                # end

                # if i < mesh.l1-1
                    # if (velocity.u.onoff[i+2] == 1.0)
                        # asum_u += AU[id, velocity.u.gIndex[i+2]] * u[velocity.u.gIndex[i+2]]
                    # end
                # end

            # end

            # pseudo_u[i] = (bU[id] - asum_u)/(AU[id,id])

        # end
    # end

    # return pseudo_u
# end

function PISO_pseudovelocity(
    AU::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,2},
    },
    AV::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,2},
    },
    velocity::CSVelocity2D,
    mesh::UnionCSMesh2D;
    T::Type{<:AbstractFloat} = Float64,
    full::Bool = true,
)
    pseudo_u = zeros(mesh.l1, mesh.m1)
    pseudo_v = zeros(mesh.l1, mesh.m1)

    ueval = phi_to_vector(velocity.u, mesh)
    veval = phi_to_vector(velocity.v, mesh)

    uiter = phi_to_vector(velocity.u, mesh; phisolution = velocity.u.iter)
    viter = phi_to_vector(velocity.v, mesh; phisolution = velocity.v.iter)

    u = ueval - uiter
    v = veval - viter

    asum_u = 0.0
    asum_v = 0.0

    for i in 1:mesh.l1
        for j in 1:mesh.m1
            if (velocity.u.onoff[i,j] == 1.0)

                id = velocity.u.gIndex[i,j]

                asum_u = 0.0
                asum_v = 0.0

                if i != 1
                    if (velocity.u.onoff[i-1,j] == 1.0)
                        asum_u += AU[id, velocity.u.gIndex[i-1,j]] * u[velocity.u.gIndex[i-1,j]]
                        asum_v += AV[id, velocity.v.gIndex[i-1,j]] * v[velocity.v.gIndex[i-1,j]]
                    end
                end

                if i != mesh.l1
                    if (velocity.u.onoff[i+1,j] == 1.0)
                        asum_u += AU[id, velocity.u.gIndex[i+1,j]] * u[velocity.u.gIndex[i+1,j]]
                        asum_v += AV[id, velocity.v.gIndex[i+1,j]] * v[velocity.v.gIndex[i+1,j]]
                    end
                end

                if j != 1
                    if (velocity.u.onoff[i,j-1] == 1.0)
                        asum_u += AU[id, velocity.u.gIndex[i,j-1]] * u[velocity.u.gIndex[i,j-1]]
                        asum_v += AV[id, velocity.v.gIndex[i,j-1]] * v[velocity.v.gIndex[i,j-1]]
                    end
                end

                if j != mesh.m1
                    if (velocity.u.onoff[i,j+1] == 1.0)
                        asum_u += AU[id, velocity.u.gIndex[i,j+1]] * u[velocity.u.gIndex[i,j+1]]
                        asum_v += AV[id, velocity.v.gIndex[i,j+1]] * v[velocity.v.gIndex[i,j+1]]
                    end
                end

                if full
                    if i > 2
                        if (velocity.u.onoff[i-2,j] == 1.0)
                            asum_u += AU[id, velocity.u.gIndex[i-2,j]] * u[velocity.u.gIndex[i-2,j]]
                            asum_v += AV[id, velocity.v.gIndex[i-2,j]] * v[velocity.v.gIndex[i-2,j]]
                        end
                    end

                    if i < mesh.l1-1
                        if (velocity.u.onoff[i+2,j] == 1.0)
                            asum_u += AU[id, velocity.u.gIndex[i+2,j]] * u[velocity.u.gIndex[i+2,j]]
                            asum_v += AV[id, velocity.v.gIndex[i+2,j]] * v[velocity.v.gIndex[i+2,j]]
                        end
                    end

                    if j > 2
                        if (velocity.u.onoff[i,j-2] == 1.0)
                            asum_u += AU[id, velocity.u.gIndex[i,j-2]] * u[velocity.u.gIndex[i,j-2]]
                            asum_v += AV[id, velocity.v.gIndex[i,j-2]] * v[velocity.v.gIndex[i,j-2]]
                        end
                    end

                    if j < mesh.m1-1
                        if (velocity.u.onoff[i,j+2] == 1.0)
                            asum_u += AU[id, velocity.u.gIndex[i,j+2]] * u[velocity.u.gIndex[i,j+2]]
                            asum_v += AV[id, velocity.v.gIndex[i,j+2]] * v[velocity.v.gIndex[i,j+2]]
                        end
                    end

                end

                pseudo_u[i,j] = -1.0 * (asum_u)/(AU[id,id])
                pseudo_v[i,j] = -1.0 * (asum_v)/(AV[id,id])

            end
        end
    end

    return pseudo_u, pseudo_v
end

# function PISO_pseudovelocity(
    # AU::Union{
    # SparseVector{<:AbstractFloat,<:Signed},
    # SparseMatrixCSC{<:AbstractFloat,<:Signed},
    # Array{<:AbstractFloat,2},
    # },
    # AV::Union{
    # SparseVector{<:AbstractFloat,<:Signed},
    # SparseMatrixCSC{<:AbstractFloat,<:Signed},
    # Array{<:AbstractFloat,2},
    # },
    # AW::Union{
    # SparseVector{<:AbstractFloat,<:Signed},
    # SparseMatrixCSC{<:AbstractFloat,<:Signed},
    # Array{<:AbstractFloat,2},
    # },
    # bU::Array{<:AbstractFloat,1},
    # bV::Array{<:AbstractFloat,1},
    # bW::Array{<:AbstractFloat,1},
    # velocity::CSVelocity3D,
    # mesh::UnionCSMesh3D;
    # relaxP = 1.0,
    # T::Type{<:AbstractFloat} = Float64,
    # full::Bool = true,
# )
    # pseudo_u = zeros(mesh.l1, mesh.m1, mesh.n1)
    # pseudo_v = zeros(mesh.l1, mesh.m1, mesh.n1)
    # pseudo_w = zeros(mesh.l1, mesh.m1, mesh.n1)

    # u = phi_to_vector(velocity.u, mesh)
    # v = phi_to_vector(velocity.v, mesh)
    # w = phi_to_vector(velocity.w, mesh)

    # asum_u = 0.0
    # asum_v = 0.0
    # asum_w = 0.0

    # for i in 1:mesh.l1
        # for j in 1:mesh.m1
            # for k in 1:mesh.n1
                # if (velocity.u.onoff[i,j,k] == 1.0)

                    # id = velocity.u.gIndex[i,j,k]

                    # asum_u = 0.0
                    # asum_v = 0.0

                    # if i != 1
                        # if (velocity.u.onoff[i-1,j] == 1.0)
                            # asum_u += AU[id, velocity.u.gIndex[i-1,j,k]] * u[velocity.u.gIndex[i-1,j,k]]
                            # asum_v += AV[id, velocity.v.gIndex[i-1,j,k]] * v[velocity.v.gIndex[i-1,j,k]]
                            # asum_w += AW[id, velocity.v.gIndex[i-1,j,k]] * w[velocity.w.gIndex[i-1,j,k]]
                        # end
                    # end

                    # if i != mesh.l1
                        # if (velocity.u.onoff[i+1,j] == 1.0)
                            # asum_u += AU[id, velocity.u.gIndex[i+1,j,k]] * u[velocity.u.gIndex[i+1,j,k]]
                            # asum_v += AV[id, velocity.v.gIndex[i+1,j,k]] * v[velocity.v.gIndex[i+1,j,k]]
                            # asum_w += AW[id, velocity.v.gIndex[i+1,j,k]] * w[velocity.w.gIndex[i+1,j,k]]
                        # end
                    # end

                    # if j != 1
                        # if (velocity.u.onoff[i,j-1,k] == 1.0)
                            # asum_u += AU[id, velocity.u.gIndex[i,j-1,k]] * u[velocity.u.gIndex[i,j-1,k]]
                            # asum_v += AV[id, velocity.v.gIndex[i,j-1,k]] * v[velocity.v.gIndex[i,j-1,k]]
                            # asum_w += AW[id, velocity.v.gIndex[i,j-1,k]] * w[velocity.w.gIndex[i,j-1,k]]
                        # end
                    # end

                    # if j != mesh.m1
                        # if (velocity.u.onoff[i,j+1,k] == 1.0)
                            # asum_u += AU[id, velocity.u.gIndex[i,j+1,k]] * u[velocity.u.gIndex[i,j+1,k]]
                            # asum_v += AV[id, velocity.v.gIndex[i,j+1,k]] * v[velocity.v.gIndex[i,j+1,k]]
                            # asum_w += AW[id, velocity.v.gIndex[i,j+1,k]] * w[velocity.w.gIndex[i,j+1,k]]
                        # end
                    # end

                    # if k != 1
                        # if (velocity.u.onoff[i,j,k-1] == 1.0)
                            # asum_u += AU[id, velocity.u.gIndex[i,j,k-1]] * u[velocity.u.gIndex[i,j,k-1]]
                            # asum_v += AV[id, velocity.v.gIndex[i,j,k-1]] * v[velocity.v.gIndex[i,j,k-1]]
                            # asum_w += AW[id, velocity.v.gIndex[i,j,k-1]] * w[velocity.w.gIndex[i,j,k-1]]
                        # end
                    # end

                    # if k != mesh.n1
                        # if (velocity.u.onoff[i,j,k+1] == 1.0)
                            # asum_u += AU[id, velocity.u.gIndex[i,j,k+1]] * u[velocity.u.gIndex[i,j,k+1]]
                            # asum_v += AV[id, velocity.v.gIndex[i,j,k+1]] * v[velocity.v.gIndex[i,j,k+1]]
                            # asum_w += AW[id, velocity.v.gIndex[i,j,k+1]] * w[velocity.w.gIndex[i,j,k+1]]
                        # end
                    # end


                    # if full
                        # if i > 2
                            # if (velocity.u.onoff[i-2,j] == 1.0)
                                # asum_u += AU[id, velocity.u.gIndex[i-2,j,k]] * u[velocity.u.gIndex[i-2,j,k]]
                                # asum_v += AV[id, velocity.v.gIndex[i-2,j,k]] * v[velocity.v.gIndex[i-2,j,k]]
                                # asum_w += AW[id, velocity.v.gIndex[i-2,j,k]] * w[velocity.w.gIndex[i-2,j,k]]
                            # end
                        # end

                        # if i < mesh.l1-1
                            # if (velocity.u.onoff[i+2,j] == 1.0)
                                # asum_u += AU[id, velocity.u.gIndex[i+2,j,k]] * u[velocity.u.gIndex[i+2,j,k]]
                                # asum_v += AV[id, velocity.v.gIndex[i+2,j,k]] * v[velocity.v.gIndex[i+2,j,k]]
                                # asum_w += AW[id, velocity.v.gIndex[i+2,j,k]] * w[velocity.w.gIndex[i+2,j,k]]
                            # end
                        # end

                        # if j > 2
                            # if (velocity.u.onoff[i,j-2,k] == 1.0)
                                # asum_u += AU[id, velocity.u.gIndex[i,j-2,k]] * u[velocity.u.gIndex[i,j-2,k]]
                                # asum_v += AV[id, velocity.v.gIndex[i,j-2,k]] * v[velocity.v.gIndex[i,j-2,k]]
                                # asum_w += AW[id, velocity.v.gIndex[i,j-2,k]] * w[velocity.w.gIndex[i,j-2,k]]
                            # end
                        # end

                        # if j < mesh.m1-1
                            # if (velocity.u.onoff[i,j+2,k] == 1.0)
                                # asum_u += AU[id, velocity.u.gIndex[i,j+2,k]] * u[velocity.u.gIndex[i,j+2,k]]
                                # asum_v += AV[id, velocity.v.gIndex[i,j+2,k]] * v[velocity.v.gIndex[i,j+2,k]]
                                # asum_w += AW[id, velocity.v.gIndex[i,j+2,k]] * w[velocity.w.gIndex[i,j+2,k]]
                            # end
                        # end

                        # if k > 2
                            # if (velocity.u.onoff[i,j,k-2] == 1.0)
                                # asum_u += AU[id, velocity.u.gIndex[i,j,k-2]] * u[velocity.u.gIndex[i,j,k-2]]
                                # asum_v += AV[id, velocity.v.gIndex[i,j,k-2]] * v[velocity.v.gIndex[i,j,k-2]]
                                # asum_w += AW[id, velocity.v.gIndex[i,j,k-2]] * w[velocity.w.gIndex[i,j,k-2]]
                            # end
                        # end

                        # if k < mesh.n1-1
                            # if (velocity.u.onoff[i,j,k+2] == 1.0)
                                # asum_u += AU[id, velocity.u.gIndex[i,j,k+2]] * u[velocity.u.gIndex[i,j,k+2]]
                                # asum_v += AV[id, velocity.v.gIndex[i,j,k+2]] * v[velocity.v.gIndex[i,j,k+2]]
                                # asum_w += AW[id, velocity.v.gIndex[i,j,k+2]] * w[velocity.w.gIndex[i,j,k+2]]
                            # end
                        # end

                    # end

                    # pseudo_u[i,j,k] = (bU[id] - asum_u)/(AU[id,id])
                    # pseudo_v[i,j,k] = (bV[id] - asum_v)/(AV[id,id])
                    # pseudo_w[i,j,k] = (bW[id] - asum_w)/(AW[id,id])

                # end
            # end
        # end
    # end

    # return pseudo_u, pseudo_v
# end
