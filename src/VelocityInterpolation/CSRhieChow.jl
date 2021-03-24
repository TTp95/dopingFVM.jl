"""

"""
function compute_RhieChow! end

function compute_RhieChow!(
    AU::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,2},
    },
    velocity::CSVelocity1D,
    mesh::UnionCSMesh1D;
    threads::Bool = false,
)
    if mesh.l1 != 1
        if threads
            Base.Threads.@threads for i in 2:mesh.l1
                if velocity.u.onoff[i-1] && velocity.u.onoff[i]
                    dx1 = 0.0
                    dx2 =  0.0
                    id1 = 0
                    id2 = 0
                    mean = 0.0
                    D1  = 0.0
                    D2  = 0.0
                    Df = 0.0
                    gradpressure = 0.0
                    a1 = 0.0
                    a2 = 0.0
                    meangradpressure = 0.0
                    pressure = 0.0

                    #parameters
                    dx1 = 0.5 * mesh.dx[i-1]
                    dx2 =  0.5 * mesh.dx[i]
                    id1 = velocity.u.gIndex[i-1]
                    id2 = velocity.u.gIndex[i]

                    #mean term
                    num = velocity.u.eval[i-1] * dx2 + velocity.u.eval[i] * dx1
                    den = dx1 + dx2
                    mean = num / den

                    #pressure term
                    D1  = mesh.vol[i-1] / AU[id1,id1]
                    D2  = mesh.vol[i] / AU[id2,id2]
                    num = D1 * dx2 + D2 * dx1
                    den = dx1 + dx2
                    Df = num / den

                    num = velocity.p.eval[i] - velocity.p.eval[i-1]
                    den = mesh.x[i] - mesh.x[i-1]
                    gradpressure = num / den

                    if (i == 2)
                        num = velocity.p.eval[i] - velocity.p.eval[i-1]
                        den = mesh.x[i] - mesh.x[i-1]
                        a1 = num / den
                    elseif !velocity.u.onoff[i-2]
                        num = velocity.p.eval[i] - velocity.p.eval[i-1]
                        den = mesh.x[i] - mesh.x[i-1]
                        a1 = num / den
                    else
                        num = velocity.p.eval[i] - velocity.p.eval[i-2]
                        den = mesh.x[i] - mesh.x[i-2]
                        a1 = num / den
                    end

                    if (i == mesh.l1)
                        num = velocity.p.eval[i] - velocity.p.eval[i-1]
                        den = mesh.x[i] - mesh.x[i-1]
                        a2 = num / den
                    elseif !velocity.u.onoff[i+1]
                        num = velocity.p.eval[i] - velocity.p.eval[i-1]
                        den = mesh.x[i] - mesh.x[i-1]
                        a2 = num / den
                    else
                        num = velocity.p.eval[i+1] - velocity.p.eval[i-1]
                        den = mesh.x[i+1] - mesh.x[i-1]
                        a2 = num / den
                    end

                    meangradpressure = 0.5 * (a1 + a2)

                    pressure = Df * (gradpressure - meangradpressure)

                    velocity.fValues.uFace[i] = mean - pressure

                elseif !velocity.u.onoff[i-1] && !velocity.u.onoff[i]
                    velocity.fValues.uFace[i] = 0.0
                end
            end

        elseif !threads
            for i in 2:mesh.l1
                if velocity.u.onoff[i-1] && velocity.u.onoff[i]
                    dx1 = 0.0
                    dx2 =  0.0
                    id1 = 0
                    id2 = 0
                    mean = 0.0
                    D1  = 0.0
                    D2  = 0.0
                    Df = 0.0
                    gradpressure = 0.0
                    a1 = 0.0
                    a2 = 0.0
                    meangradpressure = 0.0
                    pressure = 0.0

                    #parameters
                    dx1 = 0.5 * mesh.dx[i-1]
                    dx2 =  0.5 * mesh.dx[i]
                    id1 = velocity.u.gIndex[i-1]
                    id2 = velocity.u.gIndex[i]

                    #mean term
                    num = velocity.u.eval[i-1] * dx2 + velocity.u.eval[i] * dx1
                    den = dx1 + dx2
                    mean = num / den

                    #pressure term
                    D1  = mesh.vol[i-1] / AU[id1,id1]
                    D2  = mesh.vol[i] / AU[id2,id2]
                    num = D1 * dx2 + D2 * dx1
                    den = dx1 + dx2
                    Df = num / den

                    num = velocity.p.eval[i] - velocity.p.eval[i-1]
                    den = mesh.x[i] - mesh.x[i-1]
                    gradpressure = num / den

                    if (i == 2)
                        num = velocity.p.eval[i] - velocity.p.eval[i-1]
                        den = mesh.x[i] - mesh.x[i-1]
                        a1 = num / den
                    elseif !velocity.u.onoff[i-2]
                        num = velocity.p.eval[i] - velocity.p.eval[i-1]
                        den = mesh.x[i] - mesh.x[i-1]
                        a1 = num / den
                    else
                        num = velocity.p.eval[i] - velocity.p.eval[i-2]
                        den = mesh.x[i] - mesh.x[i-2]
                        a1 = num / den
                    end

                    if (i == mesh.l1)
                        num = velocity.p.eval[i] - velocity.p.eval[i-1]
                        den = mesh.x[i] - mesh.x[i-1]
                        a2 = num / den
                    elseif !velocity.u.onoff[i+1]
                        num = velocity.p.eval[i] - velocity.p.eval[i-1]
                        den = mesh.x[i] - mesh.x[i-1]
                        a2 = num / den
                    else
                        num = velocity.p.eval[i+1] - velocity.p.eval[i-1]
                        den = mesh.x[i+1] - mesh.x[i-1]
                        a2 = num / den
                    end

                    meangradpressure = 0.5 * (a1 + a2)

                    pressure = Df * (gradpressure - meangradpressure)

                    velocity.fValues.uFace[i] = mean - pressure

                elseif !velocity.u.onoff[i-1] && !velocity.u.onoff[i]
                    velocity.fValues.uFace[i] = 0.0
                end
            end

        end
    end

    return nothing
end

function compute_RhieChow!(
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
    threads::Bool = false,
)
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
                        mean = 0.0
                        D1  = 0.0
                        D2  = 0.0
                        Df = 0.0
                        gradpressure = 0.0
                        a1 = 0.0
                        a2 = 0.0
                        meangradpressure = 0.0
                        pressure = 0.0

                        #parameters
                        dx1 = 0.5 * mesh.dx[i-1]
                        dx2 =  0.5 * mesh.dx[i]
                        id1 = velocity.u.gIndex[i-1,j]
                        id2 = velocity.u.gIndex[i,j]

                        #mean term
                        num = velocity.u.eval[i-1,j] * dx2 + velocity.u.eval[i,j] * dx1
                        den = dx1 + dx2
                        mean = num / den

                        #pressure term
                        D1  = mesh.vol[i-1,j] / AU[id1,id1]
                        D2  = mesh.vol[i,j] / AU[id2,id2]
                        num = D1 * dx2 + D2 * dx1
                        den = dx1 + dx2
                        Df = num / den

                        num = velocity.p.eval[i,j] - velocity.p.eval[i-1,j]
                        den = mesh.x[i] - mesh.x[i-1]
                        gradpressure = num / den

                        if (i == 2)
                            num = velocity.p.eval[i,j] - velocity.p.eval[i-1,j]
                            den = mesh.x[i] - mesh.x[i-1]
                            a1 = num / den
                        elseif !velocity.u.onoff[i-2,j]
                            num = velocity.p.eval[i,j] - velocity.p.eval[i-1,j]
                            den = mesh.x[i] - mesh.x[i-1]
                            a1 = num / den
                        else
                            num = velocity.p.eval[i,j] - velocity.p.eval[i-2,j]
                            den = mesh.x[i] - mesh.x[i-2]
                            a1 = num / den
                        end

                        if (i == mesh.l1)
                            num = velocity.p.eval[i,j] - velocity.p.eval[i-1,j]
                            den = mesh.x[i] - mesh.x[i-1]
                            a2 = num / den
                        elseif !velocity.u.onoff[i+1,j]
                            num = velocity.p.eval[i,j] - velocity.p.eval[i-1,j]
                            den = mesh.x[i] - mesh.x[i-1]
                            a2 = num / den
                        else
                            num = velocity.p.eval[i+1,j] - velocity.p.eval[i-1,j]
                            den = mesh.x[i+1] - mesh.x[i-1]
                            a2 = num / den
                        end

                        meangradpressure = 0.5 * (a1 + a2)

                        pressure = Df * (gradpressure - meangradpressure)

                        velocity.fValues.uFace[i,j] = mean - pressure

                    elseif !velocity.u.onoff[i-1,j] && !velocity.u.onoff[i,j]
                        velocity.fValues.uFace[i,j] = 0.0
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
                        mean = 0.0
                        D1  = 0.0
                        D2  = 0.0
                        Df = 0.0
                        gradpressure = 0.0
                        a1 = 0.0
                        a2 = 0.0
                        meangradpressure = 0.0
                        pressure = 0.0

                        #parameters
                        dx1 = 0.5 * mesh.dx[i-1]
                        dx2 =  0.5 * mesh.dx[i]
                        id1 = velocity.u.gIndex[i-1,j]
                        id2 = velocity.u.gIndex[i,j]

                        #mean term
                        num = velocity.u.eval[i-1,j] * dx2 + velocity.u.eval[i,j] * dx1
                        den = dx1 + dx2
                        mean = num / den

                        #pressure term
                        D1  = mesh.vol[i-1,j] / AU[id1,id1]
                        D2  = mesh.vol[i,j] / AU[id2,id2]
                        num = D1 * dx2 + D2 * dx1
                        den = dx1 + dx2
                        Df = num / den

                        num = velocity.p.eval[i,j] - velocity.p.eval[i-1,j]
                        den = mesh.x[i] - mesh.x[i-1]
                        gradpressure = num / den

                        if (i == 2)
                            num = velocity.p.eval[i,j] - velocity.p.eval[i-1,j]
                            den = mesh.x[i] - mesh.x[i-1]
                            a1 = num / den
                        elseif !velocity.u.onoff[i-2,j]
                            num = velocity.p.eval[i,j] - velocity.p.eval[i-1,j]
                            den = mesh.x[i] - mesh.x[i-1]
                            a1 = num / den
                        else
                            num = velocity.p.eval[i,j] - velocity.p.eval[i-2,j]
                            den = mesh.x[i] - mesh.x[i-2]
                            a1 = num / den
                        end

                        if (i == mesh.l1)
                            num = velocity.p.eval[i,j] - velocity.p.eval[i-1,j]
                            den = mesh.x[i] - mesh.x[i-1]
                            a2 = num / den
                        elseif !velocity.u.onoff[i+1,j]
                            num = velocity.p.eval[i,j] - velocity.p.eval[i-1,j]
                            den = mesh.x[i] - mesh.x[i-1]
                            a2 = num / den
                        else
                            num = velocity.p.eval[i+1,j] - velocity.p.eval[i-1,j]
                            den = mesh.x[i+1] - mesh.x[i-1]
                            a2 = num / den
                        end

                        meangradpressure = 0.5 * (a1 + a2)

                        pressure = Df * (gradpressure - meangradpressure)

                        velocity.fValues.uFace[i,j] = mean - pressure

                    elseif !velocity.u.onoff[i-1,j] && !velocity.u.onoff[i,j]
                        velocity.fValues.uFace[i,j] = 0.0
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
                        mean = 0.0
                        D1  = 0.0
                        D2  = 0.0
                        Df = 0.0
                        gradpressure = 0.0
                        a1 = 0.0
                        a2 = 0.0
                        meangradpressure = 0.0
                        pressure = 0.0

                        #parameters
                        dx1 = 0.5 * mesh.dy[j-1]
                        dx2 =  0.5 * mesh.dy[j]
                        id1 = velocity.v.gIndex[i,j-1]
                        id2 = velocity.v.gIndex[i,j]

                        #mean term
                        num = velocity.v.eval[i,j-1] * dx2 + velocity.v.eval[i,j] * dx1
                        den = dx1 + dx2
                        mean = num / den

                        #pressure term
                        D1  = mesh.vol[i,j-1] / AV[id1,id1]
                        D2  = mesh.vol[i,j] / AV[id2,id2]
                        num = D1 * dx2 + D2 * dx1
                        den = dx1 + dx2
                        Df = num / den

                        num = velocity.p.eval[i,j] - velocity.p.eval[i,j-1]
                        den = mesh.y[j] - mesh.y[j-1]
                        gradpressure = num / den

                        if (j == 2)
                            num = velocity.p.eval[i,j] - velocity.p.eval[i,j-1]
                            den = mesh.y[j] - mesh.y[j-1]
                            a1 = num / den
                        elseif !velocity.v.onoff[i,j-2]
                            num = velocity.p.eval[i,j] - velocity.p.eval[i,j-1]
                            den = mesh.y[j] - mesh.y[j-1]
                            a1 = num / den
                        else
                            num = velocity.p.eval[i,j] - velocity.p.eval[i,j-2]
                            den = mesh.y[j] - mesh.y[j-2]
                            a1 = num / den
                        end

                        if (j == mesh.m1)
                            num = velocity.p.eval[i,j] - velocity.p.eval[i,j-1]
                            den = mesh.y[j] - mesh.y[j-1]
                            a2 = num / den
                        elseif !velocity.v.onoff[i,j+1]
                            num = velocity.p.eval[i,j] - velocity.p.eval[i,j-1]
                            den = mesh.y[j] - mesh.y[j-1]
                            a2 = num / den
                        else
                            num = velocity.p.eval[i,j+1] - velocity.p.eval[i,j-1]
                            den = mesh.y[j+1] - mesh.y[j-1]
                            a2 = num / den
                        end

                        meangradpressure = 0.5 * (a1 + a2)

                        pressure = Df * (gradpressure - meangradpressure)

                        velocity.fValues.vFace[i,j] = mean - pressure

                    elseif !velocity.v.onoff[i,j-1] && !velocity.v.onoff[i,j]
                        velocity.fValues.vFace[i,j] = 0.0
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
                        mean = 0.0
                        D1  = 0.0
                        D2  = 0.0
                        Df = 0.0
                        gradpressure = 0.0
                        a1 = 0.0
                        a2 = 0.0
                        meangradpressure = 0.0
                        pressure = 0.0

                        #parameters
                        dx1 = 0.5 * mesh.dy[j-1]
                        dx2 =  0.5 * mesh.dy[j]
                        id1 = velocity.v.gIndex[i,j-1]
                        id2 = velocity.v.gIndex[i,j]

                        #mean term
                        num = velocity.v.eval[i,j-1] * dx2 + velocity.v.eval[i,j] * dx1
                        den = dx1 + dx2
                        mean = num / den

                        #pressure term
                        D1  = mesh.vol[i,j-1] / AV[id1,id1]
                        D2  = mesh.vol[i,j] / AV[id2,id2]
                        num = D1 * dx2 + D2 * dx1
                        den = dx1 + dx2
                        Df = num / den

                        num = velocity.p.eval[i,j] - velocity.p.eval[i,j-1]
                        den = mesh.y[j] - mesh.y[j-1]
                        gradpressure = num / den

                        if (j == 2)
                            num = velocity.p.eval[i,j] - velocity.p.eval[i,j-1]
                            den = mesh.y[j] - mesh.y[j-1]
                            a1 = num / den
                        elseif !velocity.v.onoff[i,j-2]
                            num = velocity.p.eval[i,j] - velocity.p.eval[i,j-1]
                            den = mesh.y[j] - mesh.y[j-1]
                            a1 = num / den
                        else
                            num = velocity.p.eval[i,j] - velocity.p.eval[i,j-2]
                            den = mesh.y[j] - mesh.y[j-2]
                            a1 = num / den
                        end

                        if (j == mesh.m1)
                            num = velocity.p.eval[i,j] - velocity.p.eval[i,j-1]
                            den = mesh.y[j] - mesh.y[j-1]
                            a2 = num / den
                        elseif !velocity.v.onoff[i,j+1]
                            num = velocity.p.eval[i,j] - velocity.p.eval[i,j-1]
                            den = mesh.y[j] - mesh.y[j-1]
                            a2 = num / den
                        else
                            num = velocity.p.eval[i,j+1] - velocity.p.eval[i,j-1]
                            den = mesh.y[j+1] - mesh.y[j-1]
                            a2 = num / den
                        end

                        meangradpressure = 0.5 * (a1 + a2)

                        pressure = Df * (gradpressure - meangradpressure)

                        velocity.fValues.vFace[i,j] = mean - pressure

                    elseif !velocity.v.onoff[i,j-1] && !velocity.v.onoff[i,j]
                        velocity.fValues.vFace[i,j] = 0.0
                    end
                end
            end

        end
    end

    return nothing
end

function compute_RhieChow!(
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
    AW::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,2},
    },
    velocity::CSVelocity3D,
    mesh::UnionCSMesh3D;
    threads::Bool = false,
)
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
                            mean = 0.0
                            D1  = 0.0
                            D2  = 0.0
                            Df = 0.0
                            gradpressure = 0.0
                            a1 = 0.0
                            a2 = 0.0
                            meangradpressure = 0.0
                            pressure = 0.0

                            #parameters
                            dx1 = 0.5 * mesh.dx[i-1]
                            dx2 =  0.5 * mesh.dx[i]
                            id1 = velocity.u.gIndex[i-1,j,k]
                            id2 = velocity.u.gIndex[i,j,k]

                            #mean term
                            num = velocity.u.eval[i-1,j,k] * dx2 + velocity.u.eval[i,j,k] * dx1
                            den = dx1 + dx2
                            mean = num / den

                            #pressure term
                            D1  = mesh.vol[i-1,j,k] / AU[id1,id1]
                            D2  = mesh.vol[i,j,k] / AU[id2,id2]
                            num = D1 * dx2 + D2 * dx1
                            den = dx1 + dx2
                            Df = num / den

                            num = velocity.p.eval[i,j,k] - velocity.p.eval[i-1,j,k]
                            den = mesh.x[i] - mesh.x[i-1]
                            gradpressure = num / den

                            if (i == 2)
                                num = velocity.p.eval[i,j,k] - velocity.p.eval[i-1,j,k]
                                den = mesh.x[i] - mesh.x[i-1]
                                a1 = num / den
                            elseif !velocity.u.onoff[i-2,j,k]
                                num = velocity.p.eval[i,j,k] - velocity.p.eval[i-1,j,k]
                                den = mesh.x[i] - mesh.x[i-1]
                                a1 = num / den
                            else
                                num = velocity.p.eval[i,j,k] - velocity.p.eval[i-2,j,k]
                                den = mesh.x[i] - mesh.x[i-2]
                                a1 = num / den
                            end

                            if (i == mesh.l1)
                                num = velocity.p.eval[i,j,k] - velocity.p.eval[i-1,j,k]
                                den = mesh.x[i] - mesh.x[i-1]
                                a2 = num / den
                            elseif !velocity.u.onoff[i+1,j,k]
                                num = velocity.p.eval[i,j,k] - velocity.p.eval[i-1,j,k]
                                den = mesh.x[i] - mesh.x[i-1]
                                a2 = num / den
                            else
                                num = velocity.p.eval[i+1,j,k] - velocity.p.eval[i-1,j,k]
                                den = mesh.x[i+1] - mesh.x[i-1]
                                a2 = num / den
                            end

                            meangradpressure = 0.5 * (a1 + a2)

                            pressure = Df * (gradpressure - meangradpressure)

                            velocity.fValues.uFace[i,j,k] = mean - pressure

                        elseif !velocity.u.onoff[i-1,j,k] && !velocity.u.onoff[i,j,k]
                            velocity.fValues.uFace[i,j,k] = 0.0
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
                            mean = 0.0
                            D1  = 0.0
                            D2  = 0.0
                            Df = 0.0
                            gradpressure = 0.0
                            a1 = 0.0
                            a2 = 0.0
                            meangradpressure = 0.0
                            pressure = 0.0

                            #parameters
                            dx1 = 0.5 * mesh.dx[i-1]
                            dx2 =  0.5 * mesh.dx[i]
                            id1 = velocity.u.gIndex[i-1,j,k]
                            id2 = velocity.u.gIndex[i,j,k]

                            #mean term
                            num = velocity.u.eval[i-1,j,k] * dx2 + velocity.u.eval[i,j,k] * dx1
                            den = dx1 + dx2
                            mean = num / den

                            #pressure term
                            D1  = mesh.vol[i-1,j,k] / AU[id1,id1]
                            D2  = mesh.vol[i,j,k] / AU[id2,id2]
                            num = D1 * dx2 + D2 * dx1
                            den = dx1 + dx2
                            Df = num / den

                            num = velocity.p.eval[i,j,k] - velocity.p.eval[i-1,j,k]
                            den = mesh.x[i] - mesh.x[i-1]
                            gradpressure = num / den

                            if (i == 2)
                                num = velocity.p.eval[i,j,k] - velocity.p.eval[i-1,j,k]
                                den = mesh.x[i] - mesh.x[i-1]
                                a1 = num / den
                            elseif !velocity.u.onoff[i-2,j,k]
                                num = velocity.p.eval[i,j,k] - velocity.p.eval[i-1,j,k]
                                den = mesh.x[i] - mesh.x[i-1]
                                a1 = num / den
                            else
                                num = velocity.p.eval[i,j,k] - velocity.p.eval[i-2,j,k]
                                den = mesh.x[i] - mesh.x[i-2]
                                a1 = num / den
                            end

                            if (i == mesh.l1)
                                num = velocity.p.eval[i,j,k] - velocity.p.eval[i-1,j,k]
                                den = mesh.x[i] - mesh.x[i-1]
                                a2 = num / den
                            elseif !velocity.u.onoff[i+1,j,k]
                                num = velocity.p.eval[i,j,k] - velocity.p.eval[i-1,j,k]
                                den = mesh.x[i] - mesh.x[i-1]
                                a2 = num / den
                            else
                                num = velocity.p.eval[i+1,j,k] - velocity.p.eval[i-1,j,k]
                                den = mesh.x[i+1] - mesh.x[i-1]
                                a2 = num / den
                            end

                            meangradpressure = 0.5 * (a1 + a2)

                            pressure = Df * (gradpressure - meangradpressure)

                            velocity.fValues.uFace[i,j,k] = mean - pressure

                        elseif !velocity.u.onoff[i-1,j,k] && !velocity.u.onoff[i,j,k]
                            velocity.fValues.uFace[i,j,k] = 0.0
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
                            mean = 0.0
                            D1  = 0.0
                            D2  = 0.0
                            Df = 0.0
                            gradpressure = 0.0
                            a1 = 0.0
                            a2 = 0.0
                            meangradpressure = 0.0
                            pressure = 0.0

                            #parameters
                            dx1 = 0.5 * mesh.dy[j-1]
                            dx2 =  0.5 * mesh.dy[j]
                            id1 = velocity.v.gIndex[i,j-1,k]
                            id2 = velocity.v.gIndex[i,j,k]

                            #mean term
                            num = velocity.v.eval[i,j-1,k] * dx2 + velocity.v.eval[i,j,k] * dx1
                            den = dx1 + dx2
                            mean = num / den

                            #pressure term
                            D1  = mesh.vol[i,j-1,k] / AV[id1,id1]
                            D2  = mesh.vol[i,j,k] / AV[id2,id2]
                            num = D1 * dx2 + D2 * dx1
                            den = dx1 + dx2
                            Df = num / den

                            num = velocity.p.eval[i,j,k] - velocity.p.eval[i,j-1,k]
                            den = mesh.y[j] - mesh.y[j-1]
                            gradpressure = num / den

                            if (j == 2)
                                num = velocity.p.eval[i,j,k] - velocity.p.eval[i,j-1,k]
                                den = mesh.y[j] - mesh.y[j-1]
                                a1 = num / den
                            elseif !velocity.v.onoff[i,j-2,k]
                                num = velocity.p.eval[i,j,k] - velocity.p.eval[i,j-1,k]
                                den = mesh.y[j] - mesh.y[j-1]
                                a1 = num / den
                            else
                                num = velocity.p.eval[i,j,k] - velocity.p.eval[i,j-2,k]
                                den = mesh.y[j] - mesh.y[j-2]
                                a1 = num / den
                            end

                            if (j == mesh.m1)
                                num = velocity.p.eval[i,j,k] - velocity.p.eval[i,j-1,k]
                                den = mesh.y[j] - mesh.y[j-1]
                                a2 = num / den
                            elseif !velocity.v.onoff[i,j+1,k]
                                num = velocity.p.eval[i,j,k] - velocity.p.eval[i,j-1,k]
                                den = mesh.y[j] - mesh.y[j-1]
                                a2 = num / den
                            else
                                num = velocity.p.eval[i,j+1,k] - velocity.p.eval[i,j-1,k]
                                den = mesh.y[j+1] - mesh.y[j-1]
                                a2 = num / den
                            end

                            meangradpressure = 0.5 * (a1 + a2)

                            pressure = Df * (gradpressure - meangradpressure)

                            velocity.fValues.vFace[i,j,k] = mean - pressure

                        elseif !velocity.v.onoff[i,j-1,k] && !velocity.v.onoff[i,j,k]
                            velocity.fValues.vFace[i,j,k] = 0.0
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
                            mean = 0.0
                            D1  = 0.0
                            D2  = 0.0
                            Df = 0.0
                            gradpressure = 0.0
                            a1 = 0.0
                            a2 = 0.0
                            meangradpressure = 0.0
                            pressure = 0.0

                            #parameters
                            dx1 = 0.5 * mesh.dy[j-1]
                            dx2 =  0.5 * mesh.dy[j]
                            id1 = velocity.v.gIndex[i,j-1,k]
                            id2 = velocity.v.gIndex[i,j,k]

                            #mean term
                            num = velocity.v.eval[i,j-1,k] * dx2 + velocity.v.eval[i,j,k] * dx1
                            den = dx1 + dx2
                            mean = num / den

                            #pressure term
                            D1  = mesh.vol[i,j-1,k] / AV[id1,id1]
                            D2  = mesh.vol[i,j,k] / AV[id2,id2]
                            num = D1 * dx2 + D2 * dx1
                            den = dx1 + dx2
                            Df = num / den

                            num = velocity.p.eval[i,j,k] - velocity.p.eval[i,j-1,k]
                            den = mesh.y[j] - mesh.y[j-1]
                            gradpressure = num / den

                            if (j == 2)
                                num = velocity.p.eval[i,j,k] - velocity.p.eval[i,j-1,k]
                                den = mesh.y[j] - mesh.y[j-1]
                                a1 = num / den
                            elseif !velocity.v.onoff[i,j-2,k]
                                num = velocity.p.eval[i,j,k] - velocity.p.eval[i,j-1,k]
                                den = mesh.y[j] - mesh.y[j-1]
                                a1 = num / den
                            else
                                num = velocity.p.eval[i,j,k] - velocity.p.eval[i,j-2,k]
                                den = mesh.y[j] - mesh.y[j-2]
                                a1 = num / den
                            end

                            if (j == mesh.m1)
                                num = velocity.p.eval[i,j,k] - velocity.p.eval[i,j-1,k]
                                den = mesh.y[j] - mesh.y[j-1]
                                a2 = num / den
                            elseif !velocity.v.onoff[i,j+1,k]
                                num = velocity.p.eval[i,j,k] - velocity.p.eval[i,j-1,k]
                                den = mesh.y[j] - mesh.y[j-1]
                                a2 = num / den
                            else
                                num = velocity.p.eval[i,j+1,k] - velocity.p.eval[i,j-1,k]
                                den = mesh.y[j+1] - mesh.y[j-1]
                                a2 = num / den
                            end

                            meangradpressure = 0.5 * (a1 + a2)

                            pressure = Df * (gradpressure - meangradpressure)

                            velocity.fValues.vFace[i,j,k] = mean - pressure

                        elseif !velocity.v.onoff[i,j-1,k] && !velocity.v.onoff[i,j,k]
                            velocity.fValues.vFace[i,j,k] = 0.0
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
                            mean = 0.0
                            D1  = 0.0
                            D2  = 0.0
                            Df = 0.0
                            gradpressure = 0.0
                            a1 = 0.0
                            a2 = 0.0
                            meangradpressure = 0.0
                            pressure = 0.0

                            #parameters
                            dx1 = 0.5 * mesh.dz[k-1]
                            dx2 =  0.5 * mesh.dz[k]
                            id1 = velocity.w.gIndex[i,j,k-1]
                            id2 = velocity.w.gIndex[i,j,k]

                            #mean term
                            num = velocity.w.eval[i,j,k-1] * dx2 + velocity.w.eval[i,j,k] * dx1
                            den = dx1 + dx2
                            mean = num / den

                            #pressure term
                            D1  = mesh.vol[i,j,k-1] / AW[id1,id1]
                            D2  = mesh.vol[i,j,k] / AW[id2,id2]
                            num = D1 * dx2 + D2 * dx1
                            den = dx1 + dx2
                            Df = num / den

                            num = velocity.p.eval[i,j,k] - velocity.p.eval[i,j,k-1]
                            den = mesh.z[k] - mesh.z[k-1]
                            gradpressure = num / den

                            if (k == 2)
                                num = velocity.p.eval[i,j,k] - velocity.p.eval[i,j,k-1]
                                den = mesh.z[k] - mesh.z[k-1]
                                a1 = num / den
                            elseif !velocity.w.onoff[i,j,k-2]
                                num = velocity.p.eval[i,j,k] - velocity.p.eval[i,j,k-1]
                                den = mesh.z[k] - mesh.z[k-1]
                                a1 = num / den
                            else
                                num = velocity.p.eval[i,j,k] - velocity.p.eval[i,j,k-2]
                                den = mesh.z[k] - mesh.z[k-2]
                                a1 = num / den
                            end

                            if (k == mesh.n1)
                                num = velocity.p.eval[i,j,k] - velocity.p.eval[i,j,k-1]
                                den = mesh.z[k] - mesh.z[k-1]
                                a2 = num / den
                            elseif !velocity.w.onoff[i,j,k+1]
                                num = velocity.p.eval[i,j,k] - velocity.p.eval[i,j,k-1]
                                den = mesh.z[k] - mesh.z[k-1]
                                a2 = num / den
                            else
                                num = velocity.p.eval[i,j,k+1] - velocity.p.eval[i,j,k-1]
                                den = mesh.z[k+1] - mesh.z[k-1]
                                a2 = num / den
                            end

                            meangradpressure = 0.5 * (a1 + a2)

                            pressure = Df * (gradpressure - meangradpressure)

                            velocity.fValues.wFace[i,j,k] = mean - pressure

                        elseif !velocity.w.onoff[i,j,k-1] && !velocity.w.onoff[i,j,k]
                            velocity.fValues.wFace[i,j,k] = 0.0
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
                            mean = 0.0
                            D1  = 0.0
                            D2  = 0.0
                            Df = 0.0
                            gradpressure = 0.0
                            a1 = 0.0
                            a2 = 0.0
                            meangradpressure = 0.0
                            pressure = 0.0

                            #parameters
                            dx1 = 0.5 * mesh.dz[k-1]
                            dx2 =  0.5 * mesh.dz[k]
                            id1 = velocity.w.gIndex[i,j,k-1]
                            id2 = velocity.w.gIndex[i,j,k]

                            #mean term
                            num = velocity.w.eval[i,j,k-1] * dx2 + velocity.w.eval[i,j,k] * dx1
                            den = dx1 + dx2
                            mean = num / den

                            #pressure term
                            D1  = mesh.vol[i,j,k-1] / AW[id1,id1]
                            D2  = mesh.vol[i,j,k] / AW[id2,id2]
                            num = D1 * dx2 + D2 * dx1
                            den = dx1 + dx2
                            Df = num / den

                            num = velocity.p.eval[i,j,k] - velocity.p.eval[i,j,k-1]
                            den = mesh.z[k] - mesh.z[k-1]
                            gradpressure = num / den

                            if (k == 2)
                                num = velocity.p.eval[i,j,k] - velocity.p.eval[i,j,k-1]
                                den = mesh.z[k] - mesh.z[k-1]
                                a1 = num / den
                            elseif !velocity.w.onoff[i,j,k-2]
                                num = velocity.p.eval[i,j,k] - velocity.p.eval[i,j,k-1]
                                den = mesh.z[k] - mesh.z[k-1]
                                a1 = num / den
                            else
                                num = velocity.p.eval[i,j,k] - velocity.p.eval[i,j,k-2]
                                den = mesh.z[k] - mesh.z[k-2]
                                a1 = num / den
                            end

                            if (k == mesh.n1)
                                num = velocity.p.eval[i,j,k] - velocity.p.eval[i,j,k-1]
                                den = mesh.z[k] - mesh.z[k-1]
                                a2 = num / den
                            elseif !velocity.w.onoff[i,j,k+1]
                                num = velocity.p.eval[i,j,k] - velocity.p.eval[i,j,k-1]
                                den = mesh.z[k] - mesh.z[k-1]
                                a2 = num / den
                            else
                                num = velocity.p.eval[i,j,k+1] - velocity.p.eval[i,j,k-1]
                                den = mesh.z[k+1] - mesh.z[k-1]
                                a2 = num / den
                            end

                            meangradpressure = 0.5 * (a1 + a2)

                            pressure = Df * (gradpressure - meangradpressure)

                            velocity.fValues.wFace[i,j,k] = mean - pressure

                        elseif !velocity.w.onoff[i,j,k-1] && !velocity.w.onoff[i,j,k]
                            velocity.fValues.wFace[i,j,k] = 0.0
                        end
                    end
                end
            end


        end
    end

    return nothing
end
