"""

"""
function SIMPLEC_correction! end

function SIMPLEC_correction!(
    x::Array{<:AbstractFloat,1},
    AU::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,2},
    },
    velocity::CSVelocity1D,
    mesh::UnionCSMesh1D;
    relaxP = 1.0,
    T::Type{<:AbstractFloat} = Float64,
    threads::Bool = false,
)
    array_pc = zeros(T, mesh.l1)
    array_DU = zeros(T, mesh.l1)
    array_gradpc = zeros(T, mesh.l1)

    # vector to array
    vector_to_phi!(x, velocity.p, mesh; phisolution = array_pc, T = T, threads = threads)

    # pressure correction gradient (result = vector)
    vector_pc = pressure_phi_gradient(velocity.p, mesh; phisolution = array_pc, T = T, threads = threads)

    # vector to array (easy manipulation for next step)
    vector_to_phi!(vector_pc, velocity.p, mesh; phisolution = array_gradpc, T = T, threads = threads)

    for i in 1:mesh.l1
        if velocity.u.onoff[i]
            id = velocity.u.gIndex[i]
            array_DU[i] = 1.0/AU[id,id]
        else
            array_DU[i] = 0.0
        end
    end

    # Pressure
    velocity.p.eval = velocity.p.eval + (relaxP * array_pc)

    # Velocity
    velocity.u.eval = velocity.u.eval - (array_DU .* array_gradpc)

    # Face Velocity
    for i in 2:mesh.l1
        if velocity.u.onoff[i] && velocity.u.onoff[i-1]
            D1 = mesh.vol[i]/AU[velocity.u.gIndex[i],velocity.u.gIndex[i]]
            D2 = mesh.vol[i-1]/AU[velocity.u.gIndex[i-1],velocity.u.gIndex[i-1]]

            num = D1 * (mesh.dx[i-1] * 0.5) + D2 * (mesh.dx[i] * 0.5)
            den = 0.5 * (mesh.dx[i-1] + mesh.dx[i])
            Df = num / den

            num = array_pc[i] - array_pc[i-1]
            den = mesh.x[i] - mesh.x[i-1]
            grad_press_face = num / den

            press = Df * grad_press_face

            velocity.fValues.uFace[i] = velocity.fValues.uFace[i] - press
        end
    end

    return nothing
end

function SIMPLEC_correction!(
    x::Array{<:AbstractFloat,1},
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
    relaxP = 1.0,
    T::Type{<:AbstractFloat} = Float64,
    threads::Bool = false,
)
    array_pc = zeros(T, mesh.l1, mesh.m1)
    array_DUx = zeros(T, mesh.l1, mesh.m1)
    array_DUy = zeros(T, mesh.l1, mesh.m1)
    array_gradpcx = zeros(T, mesh.l1, mesh.m1)
    array_gradpcy = zeros(T, mesh.l1, mesh.m1)

    # vector to array
    vector_to_phi!(x, velocity.p, mesh; phisolution = array_pc, T = T, threads = threads)

    # pressure correction gradient (result = vector)
    vector_pcx, vector_pcy = pressure_phi_gradient(velocity.p, mesh; phisolution = array_pc, T = T, threads = threads)

    # vector to array (easy manipulation for next step)
    vector_to_phi!(vector_pcx, velocity.p, mesh; phisolution = array_gradpcx, T = T, threads = threads)
    vector_to_phi!(vector_pcy, velocity.p, mesh; phisolution = array_gradpcy, T = T, threads = threads)

    AUsum = sum(AU, dims = 2)
    AVsum = sum(AV, dims = 2)

    for i in 1:mesh.l1
        for j in 1:mesh.m1
            if velocity.u.onoff[i,j]
                id = velocity.u.gIndex[i,j]
                array_DUx[i,j] = 1.0/AUsum[id]
            else
                array_DUx[i,j] = 0.0
            end

            if velocity.v.onoff[i,j]
                id = velocity.v.gIndex[i,j]
                array_DUy[i,j] = 1.0/AVsum[id]
            else
                array_DUy[i,j] = 0.0
            end
        end
    end

    # Pressure
    velocity.p.eval = velocity.p.eval + (relaxP * array_pc)

    # Velocity
    velocity.u.eval = velocity.u.eval - (array_DUx .* array_gradpcx)
    velocity.v.eval = velocity.v.eval - (array_DUy .* array_gradpcy)

    # Face Velocity - U
    for i in 2:mesh.l1
        for j in 1:mesh.m1
            if velocity.u.onoff[i,j] && velocity.u.onoff[i-1,j]
                D1 = mesh.vol[i,j]/AUsum[velocity.u.gIndex[i,j]]
                D2 = mesh.vol[i-1,j]/AUsum[velocity.u.gIndex[i-1,j]]

                num = D1 * (mesh.dx[i-1] * 0.5) + D2 * (mesh.dx[i] * 0.5)
                den = 0.5 * (mesh.dx[i-1] + mesh.dx[i])
                Df = num / den

                num = array_pc[i,j] - array_pc[i-1,j]
                den = mesh.x[i] - mesh.x[i-1]
                grad_press_face = num / den

                press = Df * grad_press_face

                velocity.fValues.uFace[i,j] = velocity.fValues.uFace[i,j] - press
            end
        end
    end

    # Face Velocity - V
    for i in 1:mesh.l1
        for j in 2:mesh.m1
            if velocity.v.onoff[i,j] && velocity.v.onoff[i,j-1]
                D1 = mesh.vol[i,j]/AVsum[velocity.u.gIndex[i,j]]
                D2 = mesh.vol[i,j-1]/AVsum[velocity.u.gIndex[i,j-1]]

                num = D1 * (mesh.dy[j-1] * 0.5) + D2 * (mesh.dy[j] * 0.5)
                den = 0.5 * (mesh.dy[j-1] + mesh.dy[j])
                Df = num / den

                num = array_pc[i,j] - array_pc[i,j-1]
                den = mesh.y[j] - mesh.y[j-1]
                grad_press_face = num / den

                press = Df * grad_press_face

                velocity.fValues.vFace[i,j] = velocity.fValues.vFace[i,j] - press
            end
        end
    end


    return nothing
end

function SIMPLEC_correction!(
    x::Array{<:AbstractFloat,1},
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
    relaxP = 1.0,
    T::Type{<:AbstractFloat} = Float64,
    threads::Bool = false,
)
    array_pc = zeros(T, mesh.l1, mesh.m1, mesh.n1)
    array_DUx = zeros(T, mesh.l1, mesh.m1, mesh.n1)
    array_DUy = zeros(T, mesh.l1, mesh.m1, mesh.n1)
    array_DUz = zeros(T, mesh.l1, mesh.m1, mesh.n1)
    array_gradpcx = zeros(T, mesh.l1, mesh.m1, mesh.n1)
    array_gradpcy = zeros(T, mesh.l1, mesh.m1, mesh.n1)
    array_gradpcz = zeros(T, mesh.l1, mesh.m1, mesh.n1)

    # vector to array
    vector_to_phi!(x, velocity.p, mesh; phisolution = array_pc, T = T, threads = threads)

    # pressure correction gradient (result = vector)
    vector_pcx, vector_pcy, vector_pcz = pressure_phi_gradient(velocity.p, mesh; phisolution = array_pc, T = T, threads = threads)

    # vector to array (easy manipulation for next step)
    vector_to_phi!(vector_pcx, velocity.p, mesh; phisolution = array_gradpcx, T = T, threads = threads)
    vector_to_phi!(vector_pcy, velocity.p, mesh; phisolution = array_gradpcy, T = T, threads = threads)
    vector_to_phi!(vector_pcz, velocity.p, mesh; phisolution = array_gradpcz, T = T, threads = threads)


    for i in 1:mesh.l1
        for j in 1:mesh.m1
            for k in 1:mesh.n1
                if velocity.u.onoff[i,j,k]
                    id = velocity.u.gIndex[i,j,k]
                    array_DUx[i,j,k] = 1.0/AU[id,id]
                else
                    array_DUx[i,j,k] = 0.0
                end

                if velocity.v.onoff[i,j,k]
                    id = velocity.v.gIndex[i,j,k]
                    array_DUy[i,j,k] = 1.0/AV[id,id]
                else
                    array_DUy[i,j,k] = 0.0
                end

                if velocity.w.onoff[i,j,k]
                    id = velocity.w.gIndex[i,j,k]
                    array_DUz[i,j,k] = 1.0/AW[id,id]
                else
                    array_DUz[i,j,k] = 0.0
                end
            end
        end
    end

    # Pressure
    velocity.p.eval = velocity.p.eval + (relaxP * array_pc)

    # Velocity
    velocity.u.eval = velocity.u.eval - (array_DUx .* array_gradpcx)
    velocity.v.eval = velocity.v.eval - (array_DUy .* array_gradpcy)
    velocity.w.eval = velocity.w.eval - (array_DUz .* array_gradpcz)

    # Face Velocity - U
    for i in 2:mesh.l1
        for j in 1:mesh.m1
            for k in 1:mesh.n1
                if velocity.u.onoff[i,j,k] && velocity.u.onoff[i-1,j,k]
                    D1 = mesh.vol[i,j,k]/AU[velocity.u.gIndex[i,j,k],velocity.u.gIndex[i,j,k]]
                    D2 = mesh.vol[i-1,j,k]/AU[velocity.u.gIndex[i-1,j,k],velocity.u.gIndex[i-1,j,k]]

                    num = D1 * (mesh.dx[i-1] * 0.5) + D2 * (mesh.dx[i] * 0.5)
                    den = 0.5 * (mesh.dx[i-1] + mesh.dx[i])
                    Df = num / den

                    num = array_pc[i,j,k] - array_pc[i-1,j,k]
                    den = mesh.x[i] - mesh.x[i-1]
                    grad_press_face = num / den

                    press = Df * grad_press_face

                    velocity.fValues.uFace[i,j,k] = velocity.fValues.uFace[i,j,k] - press
                end
            end
        end
    end

    # Face Velocity - V
    for i in 1:mesh.l1
        for j in 2:mesh.m1
            for k in 1:mesh.n1
                if velocity.v.onoff[i,j,k] && velocity.v.onoff[i,j-1,k]
                    D1 = mesh.vol[i,j,k]/AV[velocity.v.gIndex[i,j,k],velocity.v.gIndex[i,j,k]]
                    D2 = mesh.vol[i,j-1,k]/AV[velocity.v.gIndex[i,j-1,k],velocity.v.gIndex[i,j-1,k]]

                    num = D1 * (mesh.dy[j-1] * 0.5) + D2 * (mesh.dy[j] * 0.5)
                    den = 0.5 * (mesh.dy[j-1] + mesh.dy[j])
                    Df = num / den

                    num = array_pc[i,j,k] - array_pc[i,j-1,k]
                    den = mesh.y[j] - mesh.y[j-1]
                    grad_press_face = num / den

                    press = Df * grad_press_face

                    velocity.fValues.vFace[i,j,k] = velocity.fValues.vFace[i,j,k] - press
                end
            end
        end
    end

    # Face Velocity - W
    for i in 1:mesh.l1
        for j in 1:mesh.m1
            for k in 2:mesh.n1
                if velocity.w.onoff[i,j,k] && velocity.w.onoff[i,j,k-1]
                    D1 = mesh.vol[i,j,k]/AW[velocity.w.gIndex[i,j,k],velocity.w.gIndex[i,j,k]]
                    D2 = mesh.vol[i,j,k-1]/AW[velocity.w.gIndex[i,j,k-1],velocity.w.gIndex[i,j,k-1]]

                    num = D1 * (mesh.dz[k-1] * 0.5) + D2 * (mesh.dz[k] * 0.5)
                    den = 0.5 * (mesh.dz[k-1] + mesh.dz[k])
                    Df = num / den

                    num = array_pc[i,j,k] - array_pc[i,j,k-1]
                    den = mesh.z[k] - mesh.z[k-1]
                    grad_press_face = num / den

                    press = Df * grad_press_face

                    velocity.fValues.wFace[i,j,k] = velocity.fValues.wFace[i,j,k] - press
                end
            end
        end
    end

    return nothing
end
