"""

"""
function compute_RhieChow_BodyForces end

function compute_RhieChow_BodyForces(
    AU::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,2},
    },
    bU::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,1},
    },
    bU_RhieChow::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,1},
    },
    velocity::CSVelocity1D,
    mesh::UnionCSMesh1D;
    T::Type{<:AbstractFloat} = Float64,
    threads::Bool = false,
)
    rhieChow_BodyForces = zeros(T, mesh.l1)

    if mesh.l1 != 1
        if threads
            Base.Threads.@threads for i in 2:mesh.l1
                if velocity.u.onoff[i-1] && velocity.u.onoff[i]
                    dx1 = 0.0
                    dx2 =  0.0
                    id1 = 0
                    id2 = 0
                    D1  = 0.0
                    D2  = 0.0
                    Df = 0.0
                    a = 0.0
                    b = 0.0
                    meanB = 0.0
                    sumBU1 = 0.0
                    sumBU2 = 0.0
                    sumBUf = 0.0


                    #parameters
                    dx1 = 0.5 * mesh.dx[i-1]
                    dx2 =  0.5 * mesh.dx[i]
                    id1 = velocity.u.gIndex[i-1]
                    id2 = velocity.u.gIndex[i]

                    D1  = mesh.vol[i-1] / AU[id1,id1]
                    D2  = mesh.vol[i] / AU[id2,id2]
                    num = D1 * dx2 + D2 * dx1
                    den = dx1 + dx2
                    Df = num / den

                    num = (bU[id1] / mesh.vol[i-1]) * dx2 + (bU[id2] / mesh.vol[i]) * dx1
                    den = dx1 + dx2
                    meanB = num / den

                    a = bU_RhieChow[id1] / mesh.vol[i-1]
                    b = bU_RhieChow[id2] / mesh.vol[i]
                    num = a * dx2 + b * dx1
                    den = dx1 + dx2
                    sumBUf = num / den

                    rhieChow_BodyForces[i] = Df * (meanB - sumBUf)
                end
            end

        elseif !threads
            for i in 2:mesh.l1
                if velocity.u.onoff[i-1] && velocity.u.onoff[i]
                    dx1 = 0.0
                    dx2 =  0.0
                    id1 = 0
                    id2 = 0
                    D1  = 0.0
                    D2  = 0.0
                    Df = 0.0
                    a = 0.0
                    b = 0.0
                    meanB = 0.0
                    sumBU1 = 0.0
                    sumBU2 = 0.0
                    sumBUf = 0.0


                    #parameters
                    dx1 = 0.5 * mesh.dx[i-1]
                    dx2 =  0.5 * mesh.dx[i]
                    id1 = velocity.u.gIndex[i-1]
                    id2 = velocity.u.gIndex[i]

                    D1  = mesh.vol[i-1] / AU[id1,id1]
                    D2  = mesh.vol[i] / AU[id2,id2]
                    num = D1 * dx2 + D2 * dx1
                    den = dx1 + dx2
                    Df = num / den

                    num = (bU[id1] / mesh.vol[i-1]) * dx2 + (bU[id2] / mesh.vol[i]) * dx1
                    den = dx1 + dx2
                    meanB = num / den

                    a = bU_RhieChow[id1] / mesh.vol[i-1]
                    b = bU_RhieChow[id2] / mesh.vol[i]
                    num = a * dx2 + b * dx1
                    den = dx1 + dx2
                    sumBUf = num / den

                    rhieChow_BodyForces[i] = Df * (meanB - sumBUf)
                end
            end
        end
    end

    return rhieChow_BodyForces
end

function compute_RhieChow_BodyForces(
    AU::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,2},
    },
    bU::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,1},
    },
    bU_RhieChow::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,1},
    },
    AV::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,2},
    },
    bV::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,1},
    },
    bV_RhieChow::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,1},
    },
    velocity::CSVelocity2D,
    mesh::UnionCSMesh2D;
    T::Type{<:AbstractFloat} = Float64,
    threads::Bool = false,
)
    u_rhieChow_BodyForces = zeros(T, mesh.l1, mesh.m1)
    v_rhieChow_BodyForces = zeros(T, mesh.l1, mesh.m1)

    #Compute u
    if mesh.l1 != 1
        if threads
            Base.Threads.@threads for i in 2:mesh.l1
                for j in 1:mesh.m1
                    if velocity.u.onoff[i-1,j] && velocity.u.onoff[i,j]
                        dx1 = 0.0
                        dx2 =  0.0
                        id1 = 0
                        id2 = 0
                        D1  = 0.0
                        D2  = 0.0
                        Df = 0.0
                        a = 0.0
                        b = 0.0
                        meanB = 0.0
                        sumBU1 = 0.0
                        sumBU2 = 0.0
                        sumBUf = 0.0

                        #parameters
                        dx1 = 0.5 * mesh.dx[i-1]
                        dx2 =  0.5 * mesh.dx[i]
                        id1 = velocity.u.gIndex[i-1,j]
                        id2 = velocity.u.gIndex[i,j]

                        D1  = mesh.vol[i-1,j] / AU[id1,id1]
                        D2  = mesh.vol[i,j] / AU[id2,id2]
                        num = D1 * dx2 + D2 * dx1
                        den = dx1 + dx2
                        Df = num / den

                        num = (bU[id1] / mesh.vol[i-1,j]) * dx2 + (bU[id2] / mesh.vol[i,j]) * dx1
                        den = dx1 + dx2
                        meanB = num / den

                        a = bU_RhieChow[id1] / mesh.vol[i-1,j]
                        b = bU_RhieChow[id2] / mesh.vol[i,j]
                        num = a * dx2 + b * dx1
                        den = dx1 + dx2
                        sumBUf = num / den

                        u_rhieChow_BodyForces[i,j] = Df * (meanB - sumBUf)
                    end
                end
            end

        elseif !threads
            for i in 2:mesh.l1
                for j in 1:mesh.m1
                    if velocity.u.onoff[i-1,j] && velocity.u.onoff[i,j]
                        dx1 = 0.0
                        dx2 =  0.0
                        id1 = 0
                        id2 = 0
                        D1  = 0.0
                        D2  = 0.0
                        Df = 0.0
                        a = 0.0
                        b = 0.0
                        meanB = 0.0
                        sumBU1 = 0.0
                        sumBU2 = 0.0
                        sumBUf = 0.0

                        #parameters
                        dx1 = 0.5 * mesh.dx[i-1]
                        dx2 =  0.5 * mesh.dx[i]
                        id1 = velocity.u.gIndex[i-1,j]
                        id2 = velocity.u.gIndex[i,j]

                        D1  = mesh.vol[i-1,j] / AU[id1,id1]
                        D2  = mesh.vol[i,j] / AU[id2,id2]
                        num = D1 * dx2 + D2 * dx1
                        den = dx1 + dx2
                        Df = num / den

                        num = (bU[id1] / mesh.vol[i-1,j]) * dx2 + (bU[id2] / mesh.vol[i,j]) * dx1
                        den = dx1 + dx2
                        meanB = num / den

                        a = bU_RhieChow[id1] / mesh.vol[i-1,j]
                        b = bU_RhieChow[id2] / mesh.vol[i,j]
                        num = a * dx2 + b * dx1
                        den = dx1 + dx2
                        sumBUf = num / den

                        u_rhieChow_BodyForces[i,j] = Df * (meanB - sumBUf)
                    end
                end
            end

        end
    end

    #Compute v
    if mesh.m1 != 1
        if threads
            Base.Threads.@threads for i in 1:mesh.l1
                for j in 2:mesh.m1
                    if velocity.v.onoff[i,j-1] && velocity.v.onoff[i,j]
                        dx1 = 0.0
                        dx2 =  0.0
                        id1 = 0
                        id2 = 0
                        D1  = 0.0
                        D2  = 0.0
                        Df = 0.0
                        a = 0.0
                        b = 0.0
                        meanB = 0.0
                        sumBU1 = 0.0
                        sumBU2 = 0.0
                        sumBUf = 0.0

                        #parameters
                        dx1 = 0.5 * mesh.dy[j-1]
                        dx2 =  0.5 * mesh.dy[j]
                        id1 = velocity.v.gIndex[i,j-1]
                        id2 = velocity.v.gIndex[i,j]

                        D1  = mesh.vol[i,j-1] / AV[id1,id1]
                        D2  = mesh.vol[i,j] / AV[id2,id2]
                        num = D1 * dx2 + D2 * dx1
                        den = dx1 + dx2
                        Df = num / den

                        num = (bV[id1] / mesh.vol[i,j-1]) * dx2 + (bV[id2] / mesh.vol[i,j]) * dx1
                        den = dx1 + dx2
                        meanB = num / den

                        a = bV_RhieChow[id1] / mesh.vol[i,j-1]
                        b = bV_RhieChow[id2] / mesh.vol[i,j]
                        num = a * dx2 + b * dx1
                        den = dx1 + dx2
                        sumBUf = num / den

                        v_rhieChow_BodyForces[i,j] = Df * (meanB - sumBUf)
                    end
                end
            end

        elseif !threads
            for i in 1:mesh.l1
                for j in 2:mesh.m1
                    if velocity.v.onoff[i,j-1] && velocity.v.onoff[i,j]
                        dx1 = 0.0
                        dx2 =  0.0
                        id1 = 0
                        id2 = 0
                        D1  = 0.0
                        D2  = 0.0
                        Df = 0.0
                        a = 0.0
                        b = 0.0
                        meanB = 0.0
                        sumBU1 = 0.0
                        sumBU2 = 0.0
                        sumBUf = 0.0

                        #parameters
                        dx1 = 0.5 * mesh.dy[j-1]
                        dx2 =  0.5 * mesh.dy[j]
                        id1 = velocity.v.gIndex[i,j-1]
                        id2 = velocity.v.gIndex[i,j]

                        D1  = mesh.vol[i,j-1] / AV[id1,id1]
                        D2  = mesh.vol[i,j] / AV[id2,id2]
                        num = D1 * dx2 + D2 * dx1
                        den = dx1 + dx2
                        Df = num / den

                        num = (bV[id1] / mesh.vol[i,j-1]) * dx2 + (bV[id2] / mesh.vol[i,j]) * dx1
                        den = dx1 + dx2
                        meanB = num / den

                        a = bV_RhieChow[id1] / mesh.vol[i,j-1]
                        b = bV_RhieChow[id2] / mesh.vol[i,j]
                        num = a * dx2 + b * dx1
                        den = dx1 + dx2
                        sumBUf = num / den

                        v_rhieChow_BodyForces[i,j] = Df * (meanB - sumBUf)
                    end
                end
            end

        end
    end

        return u_rhieChow_BodyForces, v_rhieChow_BodyForces
end

function compute_RhieChow_BodyForces(
    AU::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,2},
    },
    bU::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,1},
    },
    bU_RhieChow::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,1},
    },
    AV::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,2},
    },
    bV::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,1},
    },
    bV_RhieChow::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,1},
    },
    AW::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,2},
    },
    bW::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,1},
    },
    bW_RhieChow::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,1},
    },
    velocity::CSVelocity3D,
    mesh::UnionCSMesh3D;
    T::Type{<:AbstractFloat} = Float64,
    threads::Bool = false,
)
    u_rhieChow_BodyForces = zeros(T, mesh.l1, mesh.m1, mesh.n1)
    v_rhieChow_BodyForces = zeros(T, mesh.l1, mesh.m1, mesh.n1)
    w_rhieChow_BodyForces = zeros(T, mesh.l1, mesh.m1, mesh.n1)

    #Compute u
    if mesh.l1 != 1
        if threads
            Base.Threads.@threads for i in 2:mesh.l1
                for j in 1:mesh.m1
                    for k in 1:mesh.n1
                        if velocity.u.onoff[i-1,j,k] && velocity.u.onoff[i,j,k]
                            dx1 = 0.0
                            dx2 =  0.0
                            id1 = 0
                            id2 = 0
                            D1  = 0.0
                            D2  = 0.0
                            Df = 0.0
                            a = 0.0
                            b = 0.0
                            meanB = 0.0
                            sumBU1 = 0.0
                            sumBU2 = 0.0
                            sumBUf = 0.0

                            #parameters
                            dx1 = 0.5 * mesh.dx[i-1]
                            dx2 =  0.5 * mesh.dx[i]
                            id1 = velocity.u.gIndex[i-1,j,k]
                            id2 = velocity.u.gIndex[i,j,k]

                            D1  = mesh.vol[i-1,j,k] / AU[id1,id1]
                            D2  = mesh.vol[i,j,k] / AU[id2,id2]
                            num = D1 * dx2 + D2 * dx1
                            den = dx1 + dx2
                            Df = num / den

                            num = (bU[id1] / mesh.vol[i-1,j,k]) * dx2 + (bU[id2] / mesh.vol[i,j,k]) * dx1
                            den = dx1 + dx2
                            meanB = num / den

                            a = bU_RhieChow[id1] / mesh.vol[i-1,j,k]
                            b = bU_RhieChow[id2] / mesh.vol[i,j,k]
                            num = a * dx2 + b * dx1
                            den = dx1 + dx2
                            sumBUf = num / den

                            u_rhieChow_BodyForces[i,j,k] = Df * (meanB - sumBUf)
                        end
                    end
                end
            end

        elseif !threads
            for i in 2:mesh.l1
                for j in 1:mesh.m1
                    for k in 1:mesh.n1
                        if velocity.u.onoff[i-1,j,k] && velocity.u.onoff[i,j,k]
                            dx1 = 0.0
                            dx2 =  0.0
                            id1 = 0
                            id2 = 0
                            D1  = 0.0
                            D2  = 0.0
                            Df = 0.0
                            a = 0.0
                            b = 0.0
                            meanB = 0.0
                            sumBU1 = 0.0
                            sumBU2 = 0.0
                            sumBUf = 0.0

                            #parameters
                            dx1 = 0.5 * mesh.dx[i-1]
                            dx2 =  0.5 * mesh.dx[i]
                            id1 = velocity.u.gIndex[i-1,j,k]
                            id2 = velocity.u.gIndex[i,j,k]

                            D1  = mesh.vol[i-1,j,k] / AU[id1,id1]
                            D2  = mesh.vol[i,j,k] / AU[id2,id2]
                            num = D1 * dx2 + D2 * dx1
                            den = dx1 + dx2
                            Df = num / den

                            num = (bU[id1] / mesh.vol[i-1,j,k]) * dx2 + (bU[id2] / mesh.vol[i,j,k]) * dx1
                            den = dx1 + dx2
                            meanB = num / den

                            a = bU_RhieChow[id1] / mesh.vol[i-1,j,k]
                            b = bU_RhieChow[id2] / mesh.vol[i,j,k]
                            num = a * dx2 + b * dx1
                            den = dx1 + dx2
                            sumBUf = num / den

                            u_rhieChow_BodyForces[i,j,k] = Df * (meanB - sumBUf)
                        end
                    end
                end
            end

        end
    end

    #Compute v
    if mesh.m1 != 1
        if threads
            Base.Threads.@threads for i in 1:mesh.l1
                for j in 2:mesh.m1
                    for k in 1:mesh.n1
                        if velocity.v.onoff[i,j-1,k] && velocity.v.onoff[i,j,k]
                            dx1 = 0.0
                            dx2 =  0.0
                            id1 = 0
                            id2 = 0
                            D1  = 0.0
                            D2  = 0.0
                            Df = 0.0
                            a = 0.0
                            b = 0.0
                            meanB = 0.0
                            sumBU1 = 0.0
                            sumBU2 = 0.0
                            sumBUf = 0.0

                            #parameters
                            dx1 = 0.5 * mesh.dy[j-1]
                            dx2 =  0.5 * mesh.dy[j]
                            id1 = velocity.v.gIndex[i,j-1,k]
                            id2 = velocity.v.gIndex[i,j,k]

                            D1  = mesh.vol[i,j-1,k] / AV[id1,id1]
                            D2  = mesh.vol[i,j,k] / AV[id2,id2]
                            num = D1 * dx2 + D2 * dx1
                            den = dx1 + dx2
                            Df = num / den

                            num = (bV[id1] / mesh.vol[i,j-1,k]) * dx2 + (bV[id2] / mesh.vol[i,j,k]) * dx1
                            den = dx1 + dx2
                            meanB = num / den

                            a = bV_RhieChow[id1] / mesh.vol[i,j-1,k]
                            b = bV_RhieChow[id2] / mesh.vol[i,j,k]
                            num = a * dx2 + b * dx1
                            den = dx1 + dx2
                            sumBUf = num / den

                            v_rhieChow_BodyForces[i,j,k] = Df * (meanB - sumBUf)
                        end
                    end
                end
            end

        elseif !threads
            for i in 1:mesh.l1
                for j in 2:mesh.m1
                    for k in 1:mesh.n1
                        if velocity.v.onoff[i,j-1,k] && velocity.v.onoff[i,j,k]
                            dx1 = 0.0
                            dx2 =  0.0
                            id1 = 0
                            id2 = 0
                            D1  = 0.0
                            D2  = 0.0
                            Df = 0.0
                            a = 0.0
                            b = 0.0
                            meanB = 0.0
                            sumBU1 = 0.0
                            sumBU2 = 0.0
                            sumBUf = 0.0

                            #parameters
                            dx1 = 0.5 * mesh.dy[j-1]
                            dx2 =  0.5 * mesh.dy[j]
                            id1 = velocity.v.gIndex[i,j-1,k]
                            id2 = velocity.v.gIndex[i,j,k]

                            D1  = mesh.vol[i,j-1,k] / AV[id1,id1]
                            D2  = mesh.vol[i,j,k] / AV[id2,id2]
                            num = D1 * dx2 + D2 * dx1
                            den = dx1 + dx2
                            Df = num / den

                            num = (bV[id1] / mesh.vol[i,j-1,k]) * dx2 + (bV[id2] / mesh.vol[i,j,k]) * dx1
                            den = dx1 + dx2
                            meanB = num / den

                            a = bV_RhieChow[id1] / mesh.vol[i,j-1,k]
                            b = bV_RhieChow[id2] / mesh.vol[i,j,k]
                            num = a * dx2 + b * dx1
                            den = dx1 + dx2
                            sumBUf = num / den

                            v_rhieChow_BodyForces[i,j,k] = Df * (meanB - sumBUf)
                        end
                    end
                end
            end

        end
    end

    #Compute w
    if mesh.n1 != 1
        if threads
            Base.Threads.@threads for i in 1:mesh.l1
                for j in 1:mesh.m1
                    for k in 2:mesh.n1
                        if velocity.w.onoff[i,j,k-1] && velocity.w.onoff[i,j,k]
                            dx1 = 0.0
                            dx2 =  0.0
                            id1 = 0
                            id2 = 0
                            D1  = 0.0
                            D2  = 0.0
                            Df = 0.0
                            a = 0.0
                            b = 0.0
                            meanB = 0.0
                            sumBU1 = 0.0
                            sumBU2 = 0.0
                            sumBUf = 0.0

                            #parameters
                            dx1 = 0.5 * mesh.dz[k-1]
                            dx2 =  0.5 * mesh.dz[k]
                            id1 = velocity.w.gIndex[i,j,k-1]
                            id2 = velocity.w.gIndex[i,j,k]

                            D1  = mesh.vol[i,j,k-1] / AW[id1,id1]
                            D2  = mesh.vol[i,j,k] / AW[id2,id2]
                            num = D1 * dx2 + D2 * dx1
                            den = dx1 + dx2
                            Df = num / den

                            num = (bW[id1] / mesh.vol[i,j,k-1]) * dx2 + (bW[id2] / mesh.vol[i,j,k]) * dx1
                            den = dx1 + dx2
                            meanB = num / den

                            a = bW_RhieChow[id1] / mesh.vol[i,j,k-1]
                            b = bW_RhieChow[id2] / mesh.vol[i,j,k]
                            num = a * dx2 + b * dx1
                            den = dx1 + dx2
                            sumBUf = num / den

                            w_rhieChow_BodyForces[i,j,k] = Df * (meanB - sumBUf)
                        end
                    end
                end
            end

        elseif !threads
            for i in 1:mesh.l1
                for j in 1:mesh.m1
                    for k in 2:mesh.n1
                        if velocity.w.onoff[i,j,k-1] && velocity.w.onoff[i,j,k]
                            dx1 = 0.0
                            dx2 =  0.0
                            id1 = 0
                            id2 = 0
                            D1  = 0.0
                            D2  = 0.0
                            Df = 0.0
                            a = 0.0
                            b = 0.0
                            meanB = 0.0
                            sumBU1 = 0.0
                            sumBU2 = 0.0
                            sumBUf = 0.0

                            #parameters
                            dx1 = 0.5 * mesh.dz[k-1]
                            dx2 =  0.5 * mesh.dz[k]
                            id1 = velocity.w.gIndex[i,j,k-1]
                            id2 = velocity.w.gIndex[i,j,k]

                            D1  = mesh.vol[i,j,k-1] / AW[id1,id1]
                            D2  = mesh.vol[i,j,k] / AW[id2,id2]
                            num = D1 * dx2 + D2 * dx1
                            den = dx1 + dx2
                            Df = num / den

                            num = (bW[id1] / mesh.vol[i,j,k-1]) * dx2 + (bW[id2] / mesh.vol[i,j,k]) * dx1
                            den = dx1 + dx2
                            meanB = num / den

                            a = bW_RhieChow[id1] / mesh.vol[i,j,k-1]
                            b = bW_RhieChow[id2] / mesh.vol[i,j,k]
                            num = a * dx2 + b * dx1
                            den = dx1 + dx2
                            sumBUf = num / den

                            w_rhieChow_BodyForces[i,j,k] = Df * (meanB - sumBUf)
                        end
                    end
                end
            end

        end
    end

        return u_rhieChow_BodyForces, v_rhieChow_BodyForces, w_rhieChow_BodyForces
end
