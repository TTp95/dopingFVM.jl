"""

"""
function compute_RhieChow_Relaxation end

function compute_RhieChow_Relaxation(
    velocity::CSVelocity1D,
    mesh::UnionCSMesh1D;
    relaxU::AbstractFloat = 1.0,
    T::Type{<:AbstractFloat} = Float64,
    threads::Bool = false,
)
    rhieChow_Relax = zeros(T, mesh.l1)

    if mesh.l1 != 1
        if threads
            Base.Threads.@threads for i in 2:mesh.l1
                if velocity.u.onoff[i-1] && velocity.u.onoff[i]
                    dx1 = 0.0
                    dx2 =  0.0
                    rel = 0.0
                    velf_old = 0.0
                    meanvelf_old = 0.0
                    relaxterm = 0.0

                    #parameters
                    dx1 = 0.5 * mesh.dx[i-1]
                    dx2 =  0.5 * mesh.dx[i]

                    if relaxation
                        rel = (1.0 - relaxU)

                        velf_old = velocity.fValues.uFaceIter[i]

                        num = velocity.u.time1[i-1] * dx2 + velocity.u.time1[i] * dx1
                        den = dx1 + dx2
                        meanvelf_old = num / den

                        relaxterm = rel * (velf_old - meanvelf_old)
                    end

                    rhieChow_Relax[i] = relaxterm
                end
            end

        elseif !threads
            for i in 2:mesh.l1
                if velocity.u.onoff[i-1] && velocity.u.onoff[i]
                    dx1 = 0.0
                    dx2 =  0.0
                    rel = 0.0
                    velf_old = 0.0
                    meanvelf_old = 0.0
                    relaxterm = 0.0

                    #parameters
                    dx1 = 0.5 * mesh.dx[i-1]
                    dx2 =  0.5 * mesh.dx[i]

                    if relaxation
                        rel = (1.0 - relaxU)

                        velf_old = velocity.fValues.uFaceIter[i]

                        num = velocity.u.time1[i-1] * dx2 + velocity.u.time1[i] * dx1
                        den = dx1 + dx2
                        meanvelf_old = num / den

                        relaxterm = rel * (velf_old - meanvelf_old)
                    end

                    rhieChow_Relax[i] = relaxterm
                end
            end

        end
    end

    return rhieChow_Relax
end

function compute_RhieChow_Relaxation(
    velocity::CSVelocity2D,
    mesh::UnionCSMesh2D;
    relaxU::AbstractFloat = 1.0,
    relaxV::AbstractFloat = 1.0,
    T::Type{<:AbstractFloat} = Float64,
    threads::Bool = false,
)
    u_rhieChow_Relax = zeros(T, mesh.l1, mesh.m1)
    v_rhieChow_Relax = zeros(T, mesh.l1, mesh.m1)

    #Compute u
    if mesh.l1 != 1
        if threads
            Base.Threads.@threads for i in 2:mesh.l1
                for j in 1:mesh.m1
                    if velocity.u.onoff[i-1,j] && velocity.u.onoff[i,j]
                        dx1 = 0.0
                        dx2 =  0.0
                        rel = 0.0
                        velf_old = 0.0
                        meanvelf_old = 0.0
                        relaxterm = 0.0

                        #parameters
                        dx1 = 0.5 * mesh.dx[i-1]
                        dx2 =  0.5 * mesh.dx[i]

                        if relaxation
                            rel = (1.0 - relaxU)

                            velf_old = velocity.fValues.uFaceIter[i,j]

                            num = velocity.u.time1[i-1,j] * dx2 + velocity.u.time1[i,j] * dx1
                            den = dx1 + dx2
                            meanvelf_old = num / den

                            relaxterm = rel * (velf_old - meanvelf_old)
                        end

                        u_rhieChow_Relax[i,j] = relaxterm
                    end
                end
            end

        elseif !threads
            for i in 2:mesh.l1
                for j in 1:mesh.m1
                    if velocity.u.onoff[i-1,j] && velocity.u.onoff[i,j]
                        dx1 = 0.0
                        dx2 =  0.0
                        rel = 0.0
                        velf_old = 0.0
                        meanvelf_old = 0.0
                        relaxterm = 0.0

                        #parameters
                        dx1 = 0.5 * mesh.dx[i-1]
                        dx2 =  0.5 * mesh.dx[i]

                        if relaxation
                            rel = (1.0 - relaxU)

                            velf_old = velocity.fValues.uFaceIter[i,j]

                            num = velocity.u.time1[i-1,j] * dx2 + velocity.u.time1[i,j] * dx1
                            den = dx1 + dx2
                            meanvelf_old = num / den

                            relaxterm = rel * (velf_old - meanvelf_old)
                        end

                        u_rhieChow_Relax[i,j] = relaxterm
                    end
                end
            end

        end
    end

    #Compute v
    if mesh.m1 != 1
        if threads
            Base.Threads.@threads for i in i:mesh.l1
                for j in 2:mesh.m1
                    if velocity.v.onoff[i,j-1] && velocity.v.onoff[i,j]
                        dx1 = 0.0
                        dx2 =  0.0
                        rel = 0.0
                        velf_old = 0.0
                        meanvelf_old = 0.0
                        relaxterm = 0.0

                        #parameters
                        dx1 = 0.5 * mesh.dy[j-1]
                        dx2 =  0.5 * mesh.dy[j]

                        if relaxation
                            rel = (1.0 - relaxV)

                            velf_old = velocity.fValues.vFaceIter[i,j]

                            num = velocity.v.time1[i,j-1] * dx2 + velocity.v.time1[i,j] * dx1
                            den = dx1 + dx2
                            meanvelf_old = num / den

                            relaxterm = rel * (velf_old - meanvelf_old)
                        end

                        v_rhieChow_Relax[i,j] = relaxterm
                    end
                end
            end

        elseif !threads
            for i in i:mesh.l1
                for j in 2:mesh.m1
                    if velocity.v.onoff[i,j-1] && velocity.v.onoff[i,j]
                        dx1 = 0.0
                        dx2 =  0.0
                        rel = 0.0
                        velf_old = 0.0
                        meanvelf_old = 0.0
                        relaxterm = 0.0

                        #parameters
                        dx1 = 0.5 * mesh.dy[j-1]
                        dx2 =  0.5 * mesh.dy[j]

                        if relaxation
                            rel = (1.0 - relaxV)

                            velf_old = velocity.fValues.vFaceIter[i,j]

                            num = velocity.v.time1[i,j-1] * dx2 + velocity.v.time1[i,j] * dx1
                            den = dx1 + dx2
                            meanvelf_old = num / den

                            relaxterm = rel * (velf_old - meanvelf_old)
                        end

                        v_rhieChow_Relax[i,j] = relaxterm
                    end
                end
            end

        end
    end

    return u_rhieChow_Relax, v_rhieChow_Relax
end

function compute_RhieChow_Relaxation(
    velocity::CSVelocity3D,
    mesh::UnionCSMesh3D;
    relaxU::AbstractFloat = 1.0,
    relaxV::AbstractFloat = 1.0,
    T::Type{<:AbstractFloat} = Float64,
    threads::Bool = false,
)
    u_rhieChow_Relax = zeros(T, mesh.l1, mesh.m1, mesh.n1)
    v_rhieChow_Relax = zeros(T, mesh.l1, mesh.m1, mesh.n1)
    w_rhieChow_Relax = zeros(T, mesh.l1, mesh.m1, mesh.n1)

    #Compute u
    if mesh.l1 != 1
        if threads
            Base.Threads.@threads for i in 2:mesh.l1
                for j in 1:mesh.m1
                    for k in 1:mesh.n1
                        if velocity.u.onoff[i-1,j,k] && velocity.u.onoff[i,j,k]
                            dx1 = 0.0
                            dx2 =  0.0
                            rel = 0.0
                            velf_old = 0.0
                            meanvelf_old = 0.0
                            relaxterm = 0.0

                            #parameters
                            dx1 = 0.5 * mesh.dx[i-1]
                            dx2 =  0.5 * mesh.dx[i]

                            if relaxation
                                rel = (1.0 - relaxU)

                                velf_old = velocity.fValues.uFaceIter[i,j,k]

                                num = velocity.u.time1[i-1,j,k] * dx2 + velocity.u.time1[i,j,k] * dx1
                                den = dx1 + dx2
                                meanvelf_old = num / den

                                relaxterm = rel * (velf_old - meanvelf_old)
                            end

                            u_rhieChow_Relax[i,j,k] = relaxterm
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
                            rel = 0.0
                            velf_old = 0.0
                            meanvelf_old = 0.0
                            relaxterm = 0.0

                            #parameters
                            dx1 = 0.5 * mesh.dx[i-1]
                            dx2 =  0.5 * mesh.dx[i]

                            if relaxation
                                rel = (1.0 - relaxU)

                                velf_old = velocity.fValues.uFaceIter[i,j,k]

                                num = velocity.u.time1[i-1,j,k] * dx2 + velocity.u.time1[i,j,k] * dx1
                                den = dx1 + dx2
                                meanvelf_old = num / den

                                relaxterm = rel * (velf_old - meanvelf_old)
                            end

                            u_rhieChow_Relax[i,j,k] = relaxterm
                        end
                    end
                end
            end

        end
    end

    #Compute v
    if mesh.m1 != 1
        if threads
            Base.Threads.@threads for i in i:mesh.l1
                for j in 2:mesh.m1
                    for k in 1:mesh.n1
                        if velocity.v.onoff[i,j-1,k] && velocity.v.onoff[i,j,k]
                            dx1 = 0.0
                            dx2 =  0.0
                            rel = 0.0
                            velf_old = 0.0
                            meanvelf_old = 0.0
                            relaxterm = 0.0

                            #parameters
                            dx1 = 0.5 * mesh.dy[j-1]
                            dx2 =  0.5 * mesh.dy[j]

                            if relaxation
                                rel = (1.0 - relaxV)

                                velf_old = velocity.fValues.vFaceIter[i,j,k]

                                num = velocity.v.time1[i,j-1,k] * dx2 + velocity.v.time1[i,j,k] * dx1
                                den = dx1 + dx2
                                meanvelf_old = num / den

                                relaxterm = rel * (velf_old - meanvelf_old)
                            end

                            v_rhieChow_Relax[i,j,k] = relaxterm
                        end
                    end
                end
            end

        elseif !threads
            for i in i:mesh.l1
                for j in 2:mesh.m1
                    for k in 1:mesh.n1
                        if velocity.v.onoff[i,j-1,k] && velocity.v.onoff[i,j,k]
                            dx1 = 0.0
                            dx2 =  0.0
                            rel = 0.0
                            velf_old = 0.0
                            meanvelf_old = 0.0
                            relaxterm = 0.0

                            #parameters
                            dx1 = 0.5 * mesh.dy[j-1]
                            dx2 =  0.5 * mesh.dy[j]

                            if relaxation
                                rel = (1.0 - relaxV)

                                velf_old = velocity.fValues.vFaceIter[i,j,k]

                                num = velocity.v.time1[i,j-1,k] * dx2 + velocity.v.time1[i,j,k] * dx1
                                den = dx1 + dx2
                                meanvelf_old = num / den

                                relaxterm = rel * (velf_old - meanvelf_old)
                            end

                            v_rhieChow_Relax[i,j,k] = relaxterm
                        end
                    end
                end
            end

        end
    end

    #Compute w
    if mesh.n1 != 1
        if threads
            Base.Threads.@threads for i in i:mesh.l1
                for j in 1:mesh.m1
                    for k in 2:mesh.n1
                        if velocity.w.onoff[i,j,k-1] && velocity.w.onoff[i,j,k]
                            dx1 = 0.0
                            dx2 =  0.0
                            rel = 0.0
                            velf_old = 0.0
                            meanvelf_old = 0.0
                            relaxterm = 0.0

                            #parameters
                            dx1 = 0.5 * mesh.dz[k-1]
                            dx2 =  0.5 * mesh.dz[k]

                            if relaxation
                                rel = (1.0 - relaxW)

                                velf_old = velocity.fValues.wFaceIter[i,j,k]

                                num = velocity.w.time1[i,j,k-1] * dx2 + velocity.w.time1[i,j,k] * dx1
                                den = dx1 + dx2
                                meanvelf_old = num / den

                                relaxterm = rel * (velf_old - meanvelf_old)
                            end

                            w_rhieChow_Relax[i,j,k] = relaxterm
                        end
                    end
                end
            end

        elseif !threads
            Base.Threads.@threads for i in i:mesh.l1
                for j in 1:mesh.m1
                    for k in 2:mesh.n1
                        if velocity.w.onoff[i,j,k-1] && velocity.w.onoff[i,j,k]
                            dx1 = 0.0
                            dx2 =  0.0
                            rel = 0.0
                            velf_old = 0.0
                            meanvelf_old = 0.0
                            relaxterm = 0.0

                            #parameters
                            dx1 = 0.5 * mesh.dz[k-1]
                            dx2 =  0.5 * mesh.dz[k]

                            if relaxation
                                rel = (1.0 - relaxW)

                                velf_old = velocity.fValues.wFaceIter[i,j,k]

                                num = velocity.w.time1[i,j,k-1] * dx2 + velocity.w.time1[i,j,k] * dx1
                                den = dx1 + dx2
                                meanvelf_old = num / den

                                relaxterm = rel * (velf_old - meanvelf_old)
                            end

                            w_rhieChow_Relax[i,j,k] = relaxterm
                        end
                    end
                end
            end

        end
    end

    return u_rhieChow_Relax, v_rhieChow_Relax, w_rhieChow_Relax
end
