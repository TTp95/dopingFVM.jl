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
    pressure_value::AbstractArray = velocity.p.eval,
    mthreads::Bool = false,
)
    if !mthreads
        _compute_RhieChow!_(
            AU,
            mesh,
            velocity.u.gIndex,
            velocity.u.onoff,
            pressure_value,
            velocity.u.eval,
            velocity.fValues.uFace,
        )
    else
        _compute_RhieChow_threads!_(
            AU,
            mesh,
            velocity.u.gIndex,
            velocity.u.onoff,
            pressure_value,
            velocity.u.eval,
            velocity.fValues.uFace,
        )
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
    pressure_value::AbstractArray = velocity.p.eval,
    mthreads::Bool = false,
)
    if !mthreads
        _compute_RhieChow!_(
            AU,
            AV,
            mesh,
            velocity.u.gIndex,
            velocity.v.gIndex,
            velocity.u.onoff,
            velocity.v.onoff,
            pressure_value,
            velocity.u.eval,
            velocity.v.eval,
            velocity.fValues.uFace,
            velocity.fValues.vFace,
        )
    else
        _compute_RhieChow_threads!_(
            AU,
            AV,
            mesh,
            velocity.u.gIndex,
            velocity.v.gIndex,
            velocity.u.onoff,
            velocity.v.onoff,
            pressure_value,
            velocity.u.eval,
            velocity.v.eval,
            velocity.fValues.uFace,
            velocity.fValues.vFace,
        )
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
    pressure_value::AbstractArray = velocity.p.eval,
    mthreads::Bool = false,
)
    if !mthreads
        _compute_RhieChow!_(
            AU,
            AV,
            mesh,
            velocity.u.gIndex,
            velocity.v.gIndex,
            velocity.w.gIndex,
            velocity.u.onoff,
            velocity.v.onoff,
            velocity.w.onoff,
            pressure_value,
            velocity.u.eval,
            velocity.v.eval,
            velocity.w.eval,
            velocity.fValues.uFace,
            velocity.fValues.vFace,
            velocity.fValues.wFace,
        )
    else
        _compute_RhieChow_threads!_(
        AU,
        AV,
        mesh,
        velocity.u.gIndex,
        velocity.v.gIndex,
        velocity.w.gIndex,
        velocity.u.onoff,
        velocity.v.onoff,
        velocity.w.onoff,
        pressure_value,
        velocity.u.eval,
        velocity.v.eval,
        velocity.w.eval,
        velocity.fValues.uFace,
        velocity.fValues.vFace,
        velocity.fValues.wFace,
        )
    end

    return nothing
end

"""

"""
function _compute_RhieChow!_ end

function _compute_RhieChow!_(
    AU::AbstractArray,
    velocity::CSVelocity1D,
    mesh::UnionCSMesh1D,
    gIndex_u::AbstractArray,
    onoff_u::AbstractArray,
    value_p::AbstractArray,
    value_u::AbstractArray,
    fvel_u::AbstractArray,
)
    if mesh.l1 != 1
        for i in 2:mesh.l1
            if onoff_u[i-1] && onoff_u[i]
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
                id1 = gIndex_u[i-1]
                id2 = gIndex_u[i]

                #mean term
                num = value_u[i-1] * dx2 + value_u[i] * dx1
                den = dx1 + dx2
                mean = num / den

                #pressure term
                D1  = mesh.vol[i-1] / AU[id1,id1]
                D2  = mesh.vol[i] / AU[id2,id2]
                num = D1 * dx2 + D2 * dx1
                den = dx1 + dx2
                Df = num / den

                num = value_p[i] - value_p[i-1]
                den = mesh.x[i] - mesh.x[i-1]
                gradpressure = num / den

                if (i == 2)
                    num = value_p[i] - value_p[i-1]
                    den = mesh.x[i] - mesh.x[i-1]
                    a1 = num / den
                elseif !onoff_u[i-2]
                    num = value_p[i] - value_p[i-1]
                    den = mesh.x[i] - mesh.x[i-1]
                    a1 = num / den
                else
                    num = value_p[i] - value_p[i-2]
                    den = mesh.x[i] - mesh.x[i-2]
                    a1 = num / den
                end

                if (i == mesh.l1)
                    num = value_p[i] - value_p[i-1]
                    den = mesh.x[i] - mesh.x[i-1]
                    a2 = num / den
                elseif !onoff_u[i+1]
                    num = value_p[i] - value_p[i-1]
                    den = mesh.x[i] - mesh.x[i-1]
                    a2 = num / den
                else
                    num = value_p[i+1] - value_p[i-1]
                    den = mesh.x[i+1] - mesh.x[i-1]
                    a2 = num / den
                end

                meangradpressure = 0.5 * (a1 + a2)

                pressure = Df * (gradpressure - meangradpressure)

                fvel_u[i] = mean - pressure

            elseif !onoff_u[i-1] && !onoff_u[i]
                fvel_u[i] = 0.0
            else
                fvel_u[i] = 0.0
            end
        end
    end


    return nothing
end

function _compute_RhieChow!_(
    AU::AbstractArray,
    AV::AbstractArray,
    mesh::UnionCSMesh2D,
    gIndex_u::AbstractArray,
    gIndex_v::AbstractArray,
    onoff_u::AbstractArray,
    onoff_v::AbstractArray,
    value_p::AbstractArray,
    value_u::AbstractArray,
    value_v::AbstractArray,
    fvel_u::AbstractArray,
    fvel_v::AbstractArray,
)
    #Compute u
    if mesh.l1 != 1
        for i in 2:mesh.l1
            for j in 1:mesh.m1
                if onoff_u[i-1,j] && onoff_u[i,j]
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
                    id1 = gIndex_u[i-1,j]
                    id2 = gIndex_u[i,j]

                    #mean term
                    num = value_u[i-1,j] * dx2 + value_u[i,j] * dx1
                    den = dx1 + dx2
                    mean = num / den

                    #pressure term
                    D1  = mesh.vol[i-1,j] / AU[id1,id1]
                    D2  = mesh.vol[i,j] / AU[id2,id2]
                    num = D1 * dx2 + D2 * dx1
                    den = dx1 + dx2
                    Df = num / den

                    num = value_p[i,j] - value_p[i-1,j]
                    den = mesh.x[i] - mesh.x[i-1]
                    gradpressure = num / den

                    if (i == 2)
                        num = value_p[i,j] - value_p[i-1,j]
                        den = mesh.x[i] - mesh.x[i-1]
                        a1 = num / den
                    elseif !onoff_u[i-2,j]
                        num = value_p[i,j] - value_p[i-1,j]
                        den = mesh.x[i] - mesh.x[i-1]
                        a1 = num / den
                    else
                        num = value_p[i,j] - value_p[i-2,j]
                        den = mesh.x[i] - mesh.x[i-2]
                        a1 = num / den
                    end

                    if (i == mesh.l1)
                        num = value_p[i,j] - value_p[i-1,j]
                        den = mesh.x[i] - mesh.x[i-1]
                        a2 = num / den
                    elseif !onoff_u[i+1,j]
                        num = value_p[i,j] - value_p[i-1,j]
                        den = mesh.x[i] - mesh.x[i-1]
                        a2 = num / den
                    else
                        num = value_p[i+1,j] - value_p[i-1,j]
                        den = mesh.x[i+1] - mesh.x[i-1]
                        a2 = num / den
                    end

                    meangradpressure = 0.5 * (a1 + a2)

                    pressure = Df * (gradpressure - meangradpressure)

                    fvel_u[i,j] = mean - pressure

                elseif !onoff_u[i-1,j] && !onoff_u[i,j]
                    fvel_u[i,j] = 0.0
                else
                    fvel_u[i,j] = 0.0
                end
            end
        end
    end

    #Compute v
    if mesh.m1 != 1
        for j in 2:mesh.m1
            for i in 1:mesh.l1
                if onoff_v[i,j-1] && onoff_v[i,j]
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
                    id1 = gIndex_v[i,j-1]
                    id2 = gIndex_v[i,j]

                    #mean term
                    num = value_v[i,j-1] * dx2 + value_v[i,j] * dx1
                    den = dx1 + dx2
                    mean = num / den

                    #pressure term
                    D1  = mesh.vol[i,j-1] / AV[id1,id1]
                    D2  = mesh.vol[i,j] / AV[id2,id2]
                    num = D1 * dx2 + D2 * dx1
                    den = dx1 + dx2
                    Df = num / den

                    num = value_p[i,j] - value_p[i,j-1]
                    den = mesh.y[j] - mesh.y[j-1]
                    gradpressure = num / den

                    if (j == 2)
                        num = value_p[i,j] - value_p[i,j-1]
                        den = mesh.y[j] - mesh.y[j-1]
                        a1 = num / den
                    elseif !onoff_v[i,j-2]
                        num = value_p[i,j] - value_p[i,j-1]
                        den = mesh.y[j] - mesh.y[j-1]
                        a1 = num / den
                    else
                        num = value_p[i,j] - value_p[i,j-2]
                        den = mesh.y[j] - mesh.y[j-2]
                        a1 = num / den
                    end

                    if (j == mesh.m1)
                        num = value_p[i,j] - value_p[i,j-1]
                        den = mesh.y[j] - mesh.y[j-1]
                        a2 = num / den
                    elseif !onoff_v[i,j+1]
                        num = value_p[i,j] - value_p[i,j-1]
                        den = mesh.y[j] - mesh.y[j-1]
                        a2 = num / den
                    else
                        num = value_p[i,j+1] - value_p[i,j-1]
                        den = mesh.y[j+1] - mesh.y[j-1]
                        a2 = num / den
                    end

                    meangradpressure = 0.5 * (a1 + a2)

                    pressure = Df * (gradpressure - meangradpressure)

                    fvel_v[i,j] = mean - pressure

                elseif !onoff_v[i,j-1] && !onoff_v[i,j]
                    fvel_v[i,j] = 0.0
                else
                    fvel_v[i,j] = 0.0
                end
            end
        end
    end

    return nothing
end

function _compute_RhieChow!_(
    AU::AbstractArray,
    AV::AbstractArray,
    AW::AbstractArray,
    velocity::CSVelocity3D,
    mesh::UnionCSMesh3D,
    gIndex_u::AbstractArray,
    gIndex_v::AbstractArray,
    gIndex_w::AbstractArray,
    onoff_u::AbstractArray,
    onoff_v::AbstractArray,
    onoff_w::AbstractArray,
    value_p::AbstractArray,
    value_u::AbstractArray,
    value_v::AbstractArray,
    value_w::AbstractArray,
    fvel_u::AbstractArray,
    fvel_v::AbstractArray,
    fvel_w::AbstractArray,
)
    #Compute u
    if mesh.l1 != 1
        for k in 1:mesh.n1
            for j in 1:mesh.m1
                for i in 2:mesh.l1
                    if onoff_u[i-1,j,k] && onoff_u[i,j,k]
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
                        id1 = gIndex_u[i-1,j,k]
                        id2 = gIndex_u[i,j,k]

                        #mean term
                        num = value_u[i-1,j,k] * dx2 + value_u[i,j,k] * dx1
                        den = dx1 + dx2
                        mean = num / den

                        #pressure term
                        D1  = mesh.vol[i-1,j,k] / AU[id1,id1]
                        D2  = mesh.vol[i,j,k] / AU[id2,id2]
                        num = D1 * dx2 + D2 * dx1
                        den = dx1 + dx2
                        Df = num / den

                        num = value_p[i,j,k] - value_p[i-1,j,k]
                        den = mesh.x[i] - mesh.x[i-1]
                        gradpressure = num / den

                        if (i == 2)
                            num = value_p[i,j,k] - value_p[i-1,j,k]
                            den = mesh.x[i] - mesh.x[i-1]
                            a1 = num / den
                        elseif !onoff_u[i-2,j,k]
                            num = value_p[i,j,k] - value_p[i-1,j,k]
                            den = mesh.x[i] - mesh.x[i-1]
                            a1 = num / den
                        else
                            num = value_p[i,j,k] - value_p[i-2,j,k]
                            den = mesh.x[i] - mesh.x[i-2]
                            a1 = num / den
                        end

                        if (i == mesh.l1)
                            num = value_p[i,j,k] - value_p[i-1,j,k]
                            den = mesh.x[i] - mesh.x[i-1]
                            a2 = num / den
                        elseif !onoff_u[i+1,j,k]
                            num = value_p[i,j,k] - value_p[i-1,j,k]
                            den = mesh.x[i] - mesh.x[i-1]
                            a2 = num / den
                        else
                            num = value_p[i+1,j,k] - value_p[i-1,j,k]
                            den = mesh.x[i+1] - mesh.x[i-1]
                            a2 = num / den
                        end

                        meangradpressure = 0.5 * (a1 + a2)

                        pressure = Df * (gradpressure - meangradpressure)

                        fvel_u[i,j,k] = mean - pressure

                    elseif !onoff_u[i-1,j,k] && !onoff_u[i,j,k]
                        fvel_u[i,j,k] = 0.0
                    else
                        fvel_u[i,j,k] = 0.0
                    end
                end
            end
        end
    end

    #Compute v
    if mesh.m1 != 1
        for k in 1:mesh.n1
            for j in 2:mesh.m1
                for i in 1:mesh.l1
                    if onoff_v[i,j-1,k] && onoff_v[i,j,k]
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
                        id1 = gIndex_v[i,j-1,k]
                        id2 = gIndex_v[i,j,k]

                        #mean term
                        num = value_v[i,j-1,k] * dx2 + value_v[i,j,k] * dx1
                        den = dx1 + dx2
                        mean = num / den

                        #pressure term
                        D1  = mesh.vol[i,j-1,k] / AV[id1,id1]
                        D2  = mesh.vol[i,j,k] / AV[id2,id2]
                        num = D1 * dx2 + D2 * dx1
                        den = dx1 + dx2
                        Df = num / den

                        num = value_p[i,j,k] - value_p[i,j-1,k]
                        den = mesh.y[j] - mesh.y[j-1]
                        gradpressure = num / den

                        if (j == 2)
                            num = value_p[i,j,k] - value_p[i,j-1,k]
                            den = mesh.y[j] - mesh.y[j-1]
                            a1 = num / den
                        elseif !onoff_v[i,j-2,k]
                            num = value_p[i,j,k] - value_p[i,j-1,k]
                            den = mesh.y[j] - mesh.y[j-1]
                            a1 = num / den
                        else
                            num = value_p[i,j,k] - value_p[i,j-2,k]
                            den = mesh.y[j] - mesh.y[j-2]
                            a1 = num / den
                        end

                        if (j == mesh.m1)
                            num = value_p[i,j,k] - value_p[i,j-1,k]
                            den = mesh.y[j] - mesh.y[j-1]
                            a2 = num / den
                        elseif !onoff_v[i,j+1,k]
                            num = value_p[i,j,k] - value_p[i,j-1,k]
                            den = mesh.y[j] - mesh.y[j-1]
                            a2 = num / den
                        else
                            num = value_p[i,j+1,k] - value_p[i,j-1,k]
                            den = mesh.y[j+1] - mesh.y[j-1]
                            a2 = num / den
                        end

                        meangradpressure = 0.5 * (a1 + a2)

                        pressure = Df * (gradpressure - meangradpressure)

                        fvel_v[i,j,k] = mean - pressure

                    elseif !onoff_v[i,j-1,k] && !onoff_v[i,j,k]
                        fvel_v[i,j,k] = 0.0
                    else
                        fvel_v[i,j,k] = 0.0
                    end
                end
            end
        end
    end

    #Compute w
    if mesh.n1 != 1
        for k in 2:mesh.n1
            for j in 1:mesh.m1
                for i in 1:mesh.l1
                    if onoff_w[i,j,k-1] && onoff_w[i,j,k]
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
                        id1 = gIndex_w[i,j,k-1]
                        id2 = gIndex_w[i,j,k]

                        #mean term
                        num = value_w[i,j,k-1] * dx2 + value_w[i,j,k] * dx1
                        den = dx1 + dx2
                        mean = num / den

                        #pressure term
                        D1  = mesh.vol[i,j,k-1] / AW[id1,id1]
                        D2  = mesh.vol[i,j,k] / AW[id2,id2]
                        num = D1 * dx2 + D2 * dx1
                        den = dx1 + dx2
                        Df = num / den

                        num = value_p[i,j,k] - value_p[i,j,k-1]
                        den = mesh.z[k] - mesh.z[k-1]
                        gradpressure = num / den

                        if (k == 2)
                            num = value_p[i,j,k] - value_p[i,j,k-1]
                            den = mesh.z[k] - mesh.z[k-1]
                            a1 = num / den
                        elseif !onoff_w[i,j,k-2]
                            num = value_p[i,j,k] - value_p[i,j,k-1]
                            den = mesh.z[k] - mesh.z[k-1]
                            a1 = num / den
                        else
                            num = value_p[i,j,k] - value_p[i,j,k-2]
                            den = mesh.z[k] - mesh.z[k-2]
                            a1 = num / den
                        end

                        if (k == mesh.n1)
                            num = value_p[i,j,k] - value_p[i,j,k-1]
                            den = mesh.z[k] - mesh.z[k-1]
                            a2 = num / den
                        elseif !onoff_w[i,j,k+1]
                            num = value_p[i,j,k] - value_p[i,j,k-1]
                            den = mesh.z[k] - mesh.z[k-1]
                            a2 = num / den
                        else
                            num = value_p[i,j,k+1] - value_p[i,j,k-1]
                            den = mesh.z[k+1] - mesh.z[k-1]
                            a2 = num / den
                        end

                        meangradpressure = 0.5 * (a1 + a2)

                        pressure = Df * (gradpressure - meangradpressure)

                        fvel_w[i,j,k] = mean - pressure

                    elseif !onoff_w[i,j,k-1] && !onoff_w[i,j,k]
                        fvel_w[i,j,k] = 0.0
                    else
                        fvel_w[i,j,k] = 0.0
                    end
                end
            end
        end
    end

    return nothing
end

"""

"""
function _compute_RhieChow_threads!_ end

function _compute_RhieChow_threads!_(
    AU::AbstractArray,
    velocity::CSVelocity1D,
    mesh::UnionCSMesh1D,
    gIndex_u::AbstractArray,
    onoff_u::AbstractArray,
    value_p::AbstractArray,
    value_u::AbstractArray,
    fvel_u::AbstractArray,
)
    if mesh.l1 != 1
        Base.Threads.@threads for i in 2:mesh.l1
            if onoff_u[i-1] && onoff_u[i]
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
                id1 = gIndex_u[i-1]
                id2 = gIndex_u[i]

                #mean term
                num = value_u[i-1] * dx2 + value_u[i] * dx1
                den = dx1 + dx2
                mean = num / den

                #pressure term
                D1  = mesh.vol[i-1] / AU[id1,id1]
                D2  = mesh.vol[i] / AU[id2,id2]
                num = D1 * dx2 + D2 * dx1
                den = dx1 + dx2
                Df = num / den

                num = value_p[i] - value_p[i-1]
                den = mesh.x[i] - mesh.x[i-1]
                gradpressure = num / den

                if (i == 2)
                    num = value_p[i] - value_p[i-1]
                    den = mesh.x[i] - mesh.x[i-1]
                    a1 = num / den
                elseif !onoff_u[i-2]
                    num = value_p[i] - value_p[i-1]
                    den = mesh.x[i] - mesh.x[i-1]
                    a1 = num / den
                else
                    num = value_p[i] - value_p[i-2]
                    den = mesh.x[i] - mesh.x[i-2]
                    a1 = num / den
                end

                if (i == mesh.l1)
                    num = value_p[i] - value_p[i-1]
                    den = mesh.x[i] - mesh.x[i-1]
                    a2 = num / den
                elseif !onoff_u[i+1]
                    num = value_p[i] - value_p[i-1]
                    den = mesh.x[i] - mesh.x[i-1]
                    a2 = num / den
                else
                    num = value_p[i+1] - value_p[i-1]
                    den = mesh.x[i+1] - mesh.x[i-1]
                    a2 = num / den
                end

                meangradpressure = 0.5 * (a1 + a2)

                pressure = Df * (gradpressure - meangradpressure)

                fvel_u[i] = mean - pressure

            elseif !onoff_u[i-1] && !onoff_u[i]
                fvel_u[i] = 0.0
            else
                fvel_u[i] = 0.0
            end
        end
    end

    return nothing
end

function _compute_RhieChow_threads!_(
    AU::AbstractArray,
    AV::AbstractArray,
    mesh::UnionCSMesh2D,
    gIndex_u::AbstractArray,
    gIndex_v::AbstractArray,
    onoff_u::AbstractArray,
    onoff_v::AbstractArray,
    value_p::AbstractArray,
    value_u::AbstractArray,
    value_v::AbstractArray,
    fvel_u::AbstractArray,
    fvel_v::AbstractArray,
)
    #Compute u
    if mesh.l1 != 1
        Base.Threads.@threads for i in 2:mesh.l1
            for j in 1:mesh.m1
                if onoff_u[i-1,j] && onoff_u[i,j]
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
                    id1 = gIndex_u[i-1,j]
                    id2 = gIndex_u[i,j]

                    #mean term
                    num = value_u[i-1,j] * dx2 + value_u[i,j] * dx1
                    den = dx1 + dx2
                    mean = num / den

                    #pressure term
                    D1  = mesh.vol[i-1,j] / AU[id1,id1]
                    D2  = mesh.vol[i,j] / AU[id2,id2]
                    num = D1 * dx2 + D2 * dx1
                    den = dx1 + dx2
                    Df = num / den

                    num = value_p[i,j] - value_p[i-1,j]
                    den = mesh.x[i] - mesh.x[i-1]
                    gradpressure = num / den

                    if (i == 2)
                        num = value_p[i,j] - value_p[i-1,j]
                        den = mesh.x[i] - mesh.x[i-1]
                        a1 = num / den
                    elseif !onoff_u[i-2,j]
                        num = value_p[i,j] - value_p[i-1,j]
                        den = mesh.x[i] - mesh.x[i-1]
                        a1 = num / den
                    else
                        num = value_p[i,j] - value_p[i-2,j]
                        den = mesh.x[i] - mesh.x[i-2]
                        a1 = num / den
                    end

                    if (i == mesh.l1)
                        num = value_p[i,j] - value_p[i-1,j]
                        den = mesh.x[i] - mesh.x[i-1]
                        a2 = num / den
                    elseif !onoff_u[i+1,j]
                        num = value_p[i,j] - value_p[i-1,j]
                        den = mesh.x[i] - mesh.x[i-1]
                        a2 = num / den
                    else
                        num = value_p[i+1,j] - value_p[i-1,j]
                        den = mesh.x[i+1] - mesh.x[i-1]
                        a2 = num / den
                    end

                    meangradpressure = 0.5 * (a1 + a2)

                    pressure = Df * (gradpressure - meangradpressure)

                    fvel_u[i,j] = mean - pressure

                elseif !onoff_u[i-1,j] && !onoff_u[i,j]
                    fvel_u[i,j] = 0.0
                else
                    fvel_u[i,j] = 0.0
                end
            end
        end
    end

    #Compute v
    if mesh.m1 != 1
        Base.Threads.@threads for j in 2:mesh.m1
            for i in 1:mesh.l1
                if onoff_v[i,j-1] && onoff_v[i,j]
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
                    id1 = gIndex_v[i,j-1]
                    id2 = gIndex_v[i,j]

                    #mean term
                    num = value_v[i,j-1] * dx2 + value_v[i,j] * dx1
                    den = dx1 + dx2
                    mean = num / den

                    #pressure term
                    D1  = mesh.vol[i,j-1] / AV[id1,id1]
                    D2  = mesh.vol[i,j] / AV[id2,id2]
                    num = D1 * dx2 + D2 * dx1
                    den = dx1 + dx2
                    Df = num / den

                    num = value_p[i,j] - value_p[i,j-1]
                    den = mesh.y[j] - mesh.y[j-1]
                    gradpressure = num / den

                    if (j == 2)
                        num = value_p[i,j] - value_p[i,j-1]
                        den = mesh.y[j] - mesh.y[j-1]
                        a1 = num / den
                    elseif !onoff_v[i,j-2]
                        num = value_p[i,j] - value_p[i,j-1]
                        den = mesh.y[j] - mesh.y[j-1]
                        a1 = num / den
                    else
                        num = value_p[i,j] - value_p[i,j-2]
                        den = mesh.y[j] - mesh.y[j-2]
                        a1 = num / den
                    end

                    if (j == mesh.m1)
                        num = value_p[i,j] - value_p[i,j-1]
                        den = mesh.y[j] - mesh.y[j-1]
                        a2 = num / den
                    elseif !onoff_v[i,j+1]
                        num = value_p[i,j] - value_p[i,j-1]
                        den = mesh.y[j] - mesh.y[j-1]
                        a2 = num / den
                    else
                        num = value_p[i,j+1] - value_p[i,j-1]
                        den = mesh.y[j+1] - mesh.y[j-1]
                        a2 = num / den
                    end

                    meangradpressure = 0.5 * (a1 + a2)

                    pressure = Df * (gradpressure - meangradpressure)

                    fvel_v[i,j] = mean - pressure

                elseif !onoff_v[i,j-1] && !onoff_v[i,j]
                    fvel_v[i,j] = 0.0
                else
                    fvel_u[i,j] = 0.0
                end
            end
        end
    end

    return nothing
end


function _compute_RhieChow_threads!_(
    AU::AbstractArray,
    AV::AbstractArray,
    AW::AbstractArray,
    velocity::CSVelocity3D,
    mesh::UnionCSMesh3D,
    gIndex_u::AbstractArray,
    gIndex_v::AbstractArray,
    gIndex_w::AbstractArray,
    onoff_u::AbstractArray,
    onoff_v::AbstractArray,
    onoff_w::AbstractArray,
    value_p::AbstractArray,
    value_u::AbstractArray,
    value_v::AbstractArray,
    value_w::AbstractArray,
    fvel_u::AbstractArray,
    fvel_v::AbstractArray,
    fvel_w::AbstractArray,
)
    #Compute u
    if mesh.l1 != 1
        Base.Threads.@threads for k in 1:mesh.n1
            for j in 1:mesh.m1
                for i in 2:mesh.l1
                    if onoff_u[i-1,j,k] && onoff_u[i,j,k]
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
                        id1 = gIndex_u[i-1,j,k]
                        id2 = gIndex_u[i,j,k]

                        #mean term
                        num = value_u[i-1,j,k] * dx2 + value_u[i,j,k] * dx1
                        den = dx1 + dx2
                        mean = num / den

                        #pressure term
                        D1  = mesh.vol[i-1,j,k] / AU[id1,id1]
                        D2  = mesh.vol[i,j,k] / AU[id2,id2]
                        num = D1 * dx2 + D2 * dx1
                        den = dx1 + dx2
                        Df = num / den

                        num = value_p[i,j,k] - value_p[i-1,j,k]
                        den = mesh.x[i] - mesh.x[i-1]
                        gradpressure = num / den

                        if (i == 2)
                            num = value_p[i,j,k] - value_p[i-1,j,k]
                            den = mesh.x[i] - mesh.x[i-1]
                            a1 = num / den
                        elseif !onoff_u[i-2,j,k]
                            num = value_p[i,j,k] - value_p[i-1,j,k]
                            den = mesh.x[i] - mesh.x[i-1]
                            a1 = num / den
                        else
                            num = value_p[i,j,k] - value_p[i-2,j,k]
                            den = mesh.x[i] - mesh.x[i-2]
                            a1 = num / den
                        end

                        if (i == mesh.l1)
                            num = value_p[i,j,k] - value_p[i-1,j,k]
                            den = mesh.x[i] - mesh.x[i-1]
                            a2 = num / den
                        elseif !onoff_u[i+1,j,k]
                            num = value_p[i,j,k] - value_p[i-1,j,k]
                            den = mesh.x[i] - mesh.x[i-1]
                            a2 = num / den
                        else
                            num = value_p[i+1,j,k] - value_p[i-1,j,k]
                            den = mesh.x[i+1] - mesh.x[i-1]
                            a2 = num / den
                        end

                        meangradpressure = 0.5 * (a1 + a2)

                        pressure = Df * (gradpressure - meangradpressure)

                        fvel_u[i,j,k] = mean - pressure

                    elseif !onoff_u[i-1,j,k] && !onoff_u[i,j,k]
                        fvel_u[i,j,k] = 0.0
                    else
                        fvel_u[i,j,k] = 0.0
                    end
                end
            end
        end
    end

    #Compute v
    if mesh.m1 != 1
        Base.Threads.@threads for k in 1:mesh.n1
            for j in 2:mesh.m1
                for i in 1:mesh.l1
                    if onoff_v[i,j-1,k] && onoff_v[i,j,k]
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
                        id1 = gIndex_v[i,j-1,k]
                        id2 = gIndex_v[i,j,k]

                        #mean term
                        num = value_v[i,j-1,k] * dx2 + value_v[i,j,k] * dx1
                        den = dx1 + dx2
                        mean = num / den

                        #pressure term
                        D1  = mesh.vol[i,j-1,k] / AV[id1,id1]
                        D2  = mesh.vol[i,j,k] / AV[id2,id2]
                        num = D1 * dx2 + D2 * dx1
                        den = dx1 + dx2
                        Df = num / den

                        num = value_p[i,j,k] - value_p[i,j-1,k]
                        den = mesh.y[j] - mesh.y[j-1]
                        gradpressure = num / den

                        if (j == 2)
                            num = value_p[i,j,k] - value_p[i,j-1,k]
                            den = mesh.y[j] - mesh.y[j-1]
                            a1 = num / den
                        elseif !onoff_v[i,j-2,k]
                            num = value_p[i,j,k] - value_p[i,j-1,k]
                            den = mesh.y[j] - mesh.y[j-1]
                            a1 = num / den
                        else
                            num = value_p[i,j,k] - value_p[i,j-2,k]
                            den = mesh.y[j] - mesh.y[j-2]
                            a1 = num / den
                        end

                        if (j == mesh.m1)
                            num = value_p[i,j,k] - value_p[i,j-1,k]
                            den = mesh.y[j] - mesh.y[j-1]
                            a2 = num / den
                        elseif !onoff_v[i,j+1,k]
                            num = value_p[i,j,k] - value_p[i,j-1,k]
                            den = mesh.y[j] - mesh.y[j-1]
                            a2 = num / den
                        else
                            num = value_p[i,j+1,k] - value_p[i,j-1,k]
                            den = mesh.y[j+1] - mesh.y[j-1]
                            a2 = num / den
                        end

                        meangradpressure = 0.5 * (a1 + a2)

                        pressure = Df * (gradpressure - meangradpressure)

                        fvel_v[i,j,k] = mean - pressure

                    elseif !onoff_v[i,j-1,k] && !onoff_v[i,j,k]
                        fvel_v[i,j,k] = 0.0
                    else
                        fvel_v[i,j,k] = 0.0
                    end
                end
            end
        end
    end

    #Compute w
    if mesh.n1 != 1
        Base.Threads.@threads for k in 2:mesh.n1
            for j in 1:mesh.m1
                for i in 1:mesh.l1
                    if onoff_w[i,j,k-1] && onoff_w[i,j,k]
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
                        id1 = gIndex_w[i,j,k-1]
                        id2 = gIndex_w[i,j,k]

                        #mean term
                        num = value_w[i,j,k-1] * dx2 + value_w[i,j,k] * dx1
                        den = dx1 + dx2
                        mean = num / den

                        #pressure term
                        D1  = mesh.vol[i,j,k-1] / AW[id1,id1]
                        D2  = mesh.vol[i,j,k] / AW[id2,id2]
                        num = D1 * dx2 + D2 * dx1
                        den = dx1 + dx2
                        Df = num / den

                        num = value_p[i,j,k] - value_p[i,j,k-1]
                        den = mesh.z[k] - mesh.z[k-1]
                        gradpressure = num / den

                        if (k == 2)
                            num = value_p[i,j,k] - value_p[i,j,k-1]
                            den = mesh.z[k] - mesh.z[k-1]
                            a1 = num / den
                        elseif !onoff_w[i,j,k-2]
                            num = value_p[i,j,k] - value_p[i,j,k-1]
                            den = mesh.z[k] - mesh.z[k-1]
                            a1 = num / den
                        else
                            num = value_p[i,j,k] - value_p[i,j,k-2]
                            den = mesh.z[k] - mesh.z[k-2]
                            a1 = num / den
                        end

                        if (k == mesh.n1)
                            num = value_p[i,j,k] - value_p[i,j,k-1]
                            den = mesh.z[k] - mesh.z[k-1]
                            a2 = num / den
                        elseif !onoff_w[i,j,k+1]
                            num = value_p[i,j,k] - value_p[i,j,k-1]
                            den = mesh.z[k] - mesh.z[k-1]
                            a2 = num / den
                        else
                            num = value_p[i,j,k+1] - value_p[i,j,k-1]
                            den = mesh.z[k+1] - mesh.z[k-1]
                            a2 = num / den
                        end

                        meangradpressure = 0.5 * (a1 + a2)

                        pressure = Df * (gradpressure - meangradpressure)

                        fvel_w[i,j,k] = mean - pressure

                    elseif !onoff_w[i,j,k-1] && !onoff_w[i,j,k]
                        fvel_w[i,j,k] = 0.0
                    else
                        fvel_w[i,j,k] = 0.0
                    end
                end
            end
        end
    end

    return nothing
end
