"""

"""
function SIMPLE_correction! end

function SIMPLE_correction!(
    x::Array{<:AbstractFloat,1},
    AU::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,2},
    },
    phi::CSVelocity1D,
    mesh::UnionCSMesh1D;
    relaxP = 1.0,
    T::Type{<:AbstractFloat} = Float64,
    threads::Bool = false,
)
    array_pc = zeros(T, mesh.l1)
    array_DU = zeros(T, mesh.l1)
    array_gradpc = zeros(T, mesh.l1)

    # vector to array
    vector_to_phi!(x, phi.p, mesh; phisolution = array_pc, T = T, threads = threads)

    # pressure correction gradient (result = vector)
    vector_pc = pressure_phi_gradient(phi.p, mesh; phisolution = array_pc, T = T, threads = threads)

    # vector to array (easy manipulation for next step)
    vector_to_phi!(vector_pc, phi.p, mesh; phisolution = array_gradpc, T = T, threads = threads)

    for i in 1:mesh.l1
        if phi.u.onoff[i]
            id = phi.u.gIndex[i]
            array_DU[i] = mesh.vol[i]/AU[id,id]
        else
            array_DU[i] = 0.0
        end
    end

    # Pressure
    phi.p.eval = phi.p.eval + (relaxP * array_pc)

    # Velocity
    phi.u.eval = phi.u.eval - (array_DU .* array_gradpc)

    # Face Velocity
    for i in 2:mesh.l1
        if phi.u.onoff[i] && phi.u.onoff[i-1]
            D1 = mesh.vol[i]/AU[phi.u.gIndex[i],phi.u.gIndex[i]]
            D2 = mesh.vol[i-1]/AU[phi.u.gIndex[i-1],phi.u.gIndex[i-1]]

            num = D1 * (mesh.dx[i-1] * 0.5) + D2 * (mesh.dx[i] * 0.5)
            den = 0.5 * (mesh.dx[i-1] + mesh.dx[i])
            Df = num / den

            num = array_pc[i] - array_pc[i-1]
            den = mesh.x[i] - mesh.x[i-1]
            grad_press_face = num / den

            press = Df * grad_press_face

            phi.fValues.uFace[i] = phi.fValues.uFace[i] - press
        end
    end

    return nothing
end

function SIMPLE_correction!(
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
    phi::CSVelocity2D,
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
    vector_to_phi!(x, phi.p, mesh; phisolution = array_pc, T = T, threads = threads)

    # pressure correction gradient (result = vector)
    vector_pcx, vector_pcy = pressure_phi_gradient(phi.p, mesh; phisolution = array_pc, T = T, threads = threads)

    # vector to array (easy manipulation for next step)
    vector_to_phi!(vector_pcx, phi.p, mesh; phisolution = array_gradpcx, T = T, threads = threads)
    vector_to_phi!(vector_pcy, phi.p, mesh; phisolution = array_gradpcy, T = T, threads = threads)

    for i in 1:mesh.l1
        for j in 1:mesh.m1
            if phi.u.onoff[i,j]
                id = phi.u.gIndex[i,j]
                array_DUx[i,j] = mesh.vol[i,j]/AU[id,id]
            else
                array_DUx[i,j] = 0.0
            end

            if phi.v.onoff[i,j]
                id = phi.v.gIndex[i,j]
                array_DUy[i,j] = mesh.vol[i,j]/AV[id,id]
            else
                array_DUy[i,j] = 0.0
            end
        end
    end

    # Pressure
    phi.p.eval = phi.p.eval + (relaxP * array_pc)

    # Velocity
    phi.u.eval = phi.u.eval - (array_DUx .* array_gradpcx)
    phi.v.eval = phi.v.eval - (array_DUy .* array_gradpcy)

    # Face Velocity - U
    for i in 2:mesh.l1
        for j in 1:mesh.m1
            if phi.u.onoff[i,j] && phi.u.onoff[i-1,j]
                D1 = mesh.vol[i,j]/AU[phi.u.gIndex[i,j],phi.u.gIndex[i,j]]
                D2 = mesh.vol[i-1,j]/AU[phi.u.gIndex[i-1,j],phi.u.gIndex[i-1,j]]

                num = D1 * (mesh.dx[i-1] * 0.5) + D2 * (mesh.dx[i] * 0.5)
                den = 0.5 * (mesh.dx[i-1] + mesh.dx[i])
                Df = num / den

                num = array_pc[i,j] - array_pc[i-1,j]
                den = mesh.x[i] - mesh.x[i-1]
                grad_press_face = num / den

                press = Df * grad_press_face

                phi.fValues.uFace[i,j] = phi.fValues.uFace[i,j] - press
            end
        end
    end

    # Face Velocity - V
    for i in 1:mesh.l1
        for j in 2:mesh.m1
            if phi.v.onoff[i,j] && phi.v.onoff[i,j-1]
                D1 = mesh.vol[i,j]/AV[phi.v.gIndex[i,j],phi.v.gIndex[i,j]]
                D2 = mesh.vol[i,j-1]/AV[phi.v.gIndex[i,j-1],phi.v.gIndex[i,j-1]]

                num = D1 * (mesh.dy[j-1] * 0.5) + D2 * (mesh.dy[j] * 0.5)
                den = 0.5 * (mesh.dy[j-1] + mesh.dy[j])
                Df = num / den

                num = array_pc[i,j] - array_pc[i,j-1]
                den = mesh.y[j] - mesh.y[j-1]
                grad_press_face = num / den

                press = Df * grad_press_face

                phi.fValues.vFace[i,j] = phi.fValues.vFace[i,j] - press
            end
        end
    end

    return nothing
end

function SIMPLE_correction!(
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
    phi::CSVelocity3D,
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
    vector_to_phi!(x, phi.p, mesh; phisolution = array_pc, T = T, threads = threads)

    # pressure correction gradient (result = vector)
    vector_pcx, vector_pcy, vector_pcz = pressure_phi_gradient(phi.p, mesh; phisolution = array_pc, T = T, threads = threads)

    # vector to array (easy manipulation for next step)
    vector_to_phi!(vector_pcx, phi.p, mesh; phisolution = array_gradpcx, T = T, threads = threads)
    vector_to_phi!(vector_pcy, phi.p, mesh; phisolution = array_gradpcy, T = T, threads = threads)
    vector_to_phi!(vector_pcz, phi.p, mesh; phisolution = array_gradpcz, T = T, threads = threads)


    for i in 1:mesh.l1
        for j in 1:mesh.m1
            for k in 1:mesh.n1
                if phi.u.onoff[i,j,k]
                    id = phi.u.gIndex[i,j,k]
                    array_DUx[i,j,k] = mesh.vol[i,j,k]/AU[id,id]
                else
                    array_DUx[i,j,k] = 0.0
                end

                if phi.v.onoff[i,j,k]
                    id = phi.v.gIndex[i,j,k]
                    array_DUy[i,j,k] = mesh.vol[i,j,k]/AV[id,id]
                else
                    array_DUy[i,j,k] = 0.0
                end

                if phi.w.onoff[i,j,k]
                    id = phi.w.gIndex[i,j,k]
                    array_DUz[i,j,k] = mesh.vol[i,j,k]/AW[id,id]
                else
                    array_DUz[i,j,k] = 0.0
                end
            end
        end
    end

    # Pressure
    phi.p.eval = phi.p.eval + (relaxP * array_pc)

    # Velocity
    phi.u.eval = phi.u.eval - (array_DUx .* array_gradpcx)
    phi.v.eval = phi.v.eval - (array_DUy .* array_gradpcy)
    phi.w.eval = phi.w.eval - (array_DUz .* array_gradpcz)

    # Face Velocity - U
    for i in 2:mesh.l1
        for j in 1:mesh.m1
            for k in 1:mesh.n1
                if phi.u.onoff[i,j,k] && phi.u.onoff[i-1,j,k]
                    D1 = mesh.vol[i,j,k]/AU[phi.u.gIndex[i,j,k],phi.u.gIndex[i,j,k]]
                    D2 = mesh.vol[i-1,j,k]/AU[phi.u.gIndex[i-1,j,k],phi.u.gIndex[i-1,j,k]]

                    num = D1 * (mesh.dx[i-1] * 0.5) + D2 * (mesh.dx[i] * 0.5)
                    den = 0.5 * (mesh.dx[i-1] + mesh.dx[i])
                    Df = num / den

                    num = array_pc[i,j,k] - array_pc[i-1,j,k]
                    den = mesh.x[i] - mesh.x[i-1]
                    grad_press_face = num / den

                    press = Df * grad_press_face

                    phi.fValues.uFace[i,j,k] = phi.fValues.uFace[i,j,k] - press
                end
            end
        end
    end

    # Face Velocity - V
    for i in 1:mesh.l1
        for j in 2:mesh.m1
            for k in 1:mesh.n1
                if phi.v.onoff[i,j,k] && phi.v.onoff[i,j-1,k]
                    D1 = mesh.vol[i,j,k]/AV[phi.v.gIndex[i,j,k],phi.v.gIndex[i,j,k]]
                    D2 = mesh.vol[i,j-1,k]/AV[phi.v.gIndex[i,j-1,k],phi.v.gIndex[i,j-1,k]]

                    num = D1 * (mesh.dy[j-1] * 0.5) + D2 * (mesh.dy[j] * 0.5)
                    den = 0.5 * (mesh.dy[j-1] + mesh.dy[j])
                    Df = num / den

                    num = array_pc[i,j,k] - array_pc[i,j-1,k]
                    den = mesh.y[j] - mesh.y[j-1]
                    grad_press_face = num / den

                    press = Df * grad_press_face

                    phi.fValues.vFace[i,j,k] = phi.fValues.vFace[i,j,k] - press
                end
            end
        end
    end

    # Face Velocity - W
    for i in 1:mesh.l1
        for j in 1:mesh.m1
            for k in 2:mesh.n1
                if phi.w.onoff[i,j,k] && phi.w.onoff[i,j,k-1]
                    D1 = mesh.vol[i,j,k]/AW[phi.w.gIndex[i,j,k],phi.w.gIndex[i,j,k]]
                    D2 = mesh.vol[i,j,k-1]/AW[phi.w.gIndex[i,j,k-1],phi.w.gIndex[i,j,k-1]]

                    num = D1 * (mesh.dz[k-1] * 0.5) + D2 * (mesh.dz[k] * 0.5)
                    den = 0.5 * (mesh.dz[k-1] + mesh.dz[k])
                    Df = num / den

                    num = array_pc[i,j,k] - array_pc[i,j,k-1]
                    den = mesh.z[k] - mesh.z[k-1]
                    grad_press_face = num / den

                    press = Df * grad_press_face

                    phi.fValues.wFace[i,j,k] = phi.fValues.wFace[i,j,k] - press
                end
            end
        end
    end

    return nothing
end
