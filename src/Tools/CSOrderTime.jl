"""

"""
function order_timeCS! end

function order_timeCS!(
    phi::CSPhi1D,
    mesh::UnionCSMesh1D;
    threads = false,
)

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            @inbounds phi.time3[i] = phi.time2[i]
            @inbounds phi.time2[i] = phi.time1[i]
            @inbounds phi.time1[i] = phi.eval[i]
        end
    elseif !threads
        for i in 1:mesh.l1
            @inbounds phi.time3[i] = phi.time2[i]
            @inbounds phi.time2[i] = phi.time1[i]
            @inbounds phi.time1[i] = phi.eval[i]
        end
    end

    return nothing
end

function order_timeCS!(
    phi::CSPhi2D,
    mesh::UnionCSMesh2D;
    threads = false,
)

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            for j in 1:mesh.m1
                @inbounds phi.time3[i,j] = phi.time2[i,j]
                @inbounds phi.time2[i,j] = phi.time1[i,j]
                @inbounds phi.time1[i,j] = phi.eval[i,j]
            end
        end
    elseif !threads
        for i in 1:mesh.l1
            for j in 1:mesh.m1
                @inbounds phi.time3[i,j] = phi.time2[i,j]
                @inbounds phi.time2[i,j] = phi.time1[i,j]
                @inbounds phi.time1[i,j] = phi.eval[i,j]
            end
        end
    end

    return nothing
end

function order_timeCS!(
    phi::CSPhi3D,
    mesh::UnionCSMesh3D;
    threads = false,
)

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            for j in 1:mesh.m1
                for k in 1:mesh.n1
                    @inbounds phi.time3[i,j,k] = phi.time2[i,j,k]
                    @inbounds phi.time2[i,j,k] = phi.time1[i,j,k]
                    @inbounds phi.time1[i,j,k] = phi.eval[i,j,k]
                end
            end
        end
    elseif !threads
        for i in 1:mesh.l1
            for j in 1:mesh.m1
                for k in 1:mesh.n1
                    @inbounds phi.time3[i,j,k] = phi.time2[i,j,k]
                    @inbounds phi.time2[i,j,k] = phi.time1[i,j,k]
                    @inbounds phi.time1[i,j,k] = phi.eval[i,j,k]
                end
            end
        end
    end

    return nothing
end

function order_timeCS!(
    phi::CSVelocity1D,
    mesh::UnionCSMesh1D;
    threads = false,
)
    order_timeCS!(phi.u, mesh; threads=threads)
    order_timeCS!(phi.p, mesh; threads=threads)

    if threads
        Base.Threads.@threads for i in 1:(mesh.l1 + 1)
            @inbounds phi.fValues.uFaceTime3[i] = phi.fValues.uFaceTime2[i]
            @inbounds phi.fValues.uFaceTime2[i] = phi.fValues.uFaceTime1[i]
            @inbounds phi.fValues.uFaceTime1[i] = phi.fValues.uFace[i]
        end
    elseif !threads
        for i in 1:(mesh.l1 + 1)
            @inbounds phi.fValues.uFaceTime3[i] = phi.fValues.uFaceTime2[i]
            @inbounds phi.fValues.uFaceTime2[i] = phi.fValues.uFaceTime1[i]
            @inbounds phi.fValues.uFaceTime1[i] = phi.fValues.uFace[i]
        end
    end

    return nothing
end

function order_timeCS!(
    phi::CSVelocity2D,
    mesh::UnionCSMesh2D;
    threads = false,
)
    order_timeCS!(phi.u, mesh; threads=threads)
    order_timeCS!(phi.v, mesh; threads=threads)
    order_timeCS!(phi.p, mesh; threads=threads)

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            for j in 1:mesh.m1
                if (i == mesh.m1)
                    @inbounds phi.fValues.uFaceTime3[i,j] = phi.fValues.uFaceTime2[i,j]
                    @inbounds phi.fValues.uFaceTime2[i,j] = phi.fValues.uFaceTime1[i,j]
                    @inbounds phi.fValues.uFaceTime1[i,j] = phi.fValues.uFace[i,j]
                    @inbounds phi.fValues.uFaceTime3[i+1,j] = phi.fValues.uFaceTime2[i+1,j]
                    @inbounds phi.fValues.uFaceTime2[i+1,j] = phi.fValues.uFaceTime1[i+1,j]
                    @inbounds phi.fValues.uFaceTime1[i+1,j] = phi.fValues.uFace[i+1,j]
                else
                    @inbounds phi.fValues.uFaceTime3[i,j] = phi.fValues.uFaceTime2[i,j]
                    @inbounds phi.fValues.uFaceTime2[i,j] = phi.fValues.uFaceTime1[i,j]
                    @inbounds phi.fValues.uFaceTime1[i,j] = phi.fValues.uFace[i,j]
                end
                if (i == mesh.l1)
                    @inbounds phi.fValues.vFaceTime3[i,j] = phi.fValues.vFaceTime2[i,j]
                    @inbounds phi.fValues.vFaceTime2[i,j] = phi.fValues.vFaceTime1[i,j]
                    @inbounds phi.fValues.vFaceTime1[i,j] = phi.fValues.vFace[i,j]
                    @inbounds phi.fValues.vFaceTime3[i,j+1] = phi.fValues.vFaceTime2[i,j+1]
                    @inbounds phi.fValues.vFaceTime2[i,j+1] = phi.fValues.vFaceTime1[i,j+1]
                    @inbounds phi.fValues.vFaceTime1[i,j+1] = phi.fValues.vFace[i,j+1]
                else
                    @inbounds phi.fValues.vFaceTime3[i,j] = phi.fValues.vFaceTime2[i,j]
                    @inbounds phi.fValues.vFaceTime2[i,j] = phi.fValues.vFaceTime1[i,j]
                    @inbounds phi.fValues.vFaceTime1[i,j] = phi.fValues.vFace[i,j]
                end
            end
        end
    elseif !threads
        for i in 1:mesh.l1
            for j in 1:mesh.m1
                if (i == mesh.m1)
                    @inbounds phi.fValues.uFaceTime3[i,j] = phi.fValues.uFaceTime2[i,j]
                    @inbounds phi.fValues.uFaceTime2[i,j] = phi.fValues.uFaceTime1[i,j]
                    @inbounds phi.fValues.uFaceTime1[i,j] = phi.fValues.uFace[i,j]
                    @inbounds phi.fValues.uFaceTime3[i+1,j] = phi.fValues.uFaceTime2[i+1,j]
                    @inbounds phi.fValues.uFaceTime2[i+1,j] = phi.fValues.uFaceTime1[i+1,j]
                    @inbounds phi.fValues.uFaceTime1[i+1,j] = phi.fValues.uFace[i+1,j]
                else
                    @inbounds phi.fValues.uFaceTime3[i,j] = phi.fValues.uFaceTime2[i,j]
                    @inbounds phi.fValues.uFaceTime2[i,j] = phi.fValues.uFaceTime1[i,j]
                    @inbounds phi.fValues.uFaceTime1[i,j] = phi.fValues.uFace[i,j]
                end
                if (i == mesh.l1)
                    @inbounds phi.fValues.vFaceTime3[i,j] = phi.fValues.vFaceTime2[i,j]
                    @inbounds phi.fValues.vFaceTime2[i,j] = phi.fValues.vFaceTime1[i,j]
                    @inbounds phi.fValues.vFaceTime1[i,j] = phi.fValues.vFace[i,j]
                    @inbounds phi.fValues.vFaceTime3[i,j+1] = phi.fValues.vFaceTime2[i,j+1]
                    @inbounds phi.fValues.vFaceTime2[i,j+1] = phi.fValues.vFaceTime1[i,j+1]
                    @inbounds phi.fValues.vFaceTime1[i,j+1] = phi.fValues.vFace[i,j+1]
                else
                    @inbounds phi.fValues.vFaceTime3[i,j] = phi.fValues.vFaceTime2[i,j]
                    @inbounds phi.fValues.vFaceTime2[i,j] = phi.fValues.vFaceTime1[i,j]
                    @inbounds phi.fValues.vFaceTime1[i,j] = phi.fValues.vFace[i,j]
                end
            end
        end
    end

    return nothing
end

function order_timeCS!(
    phi::CSVelocity3D,
    mesh::UnionCSMesh3D;
    threads = false,
)
    order_timeCS!(phi.u, mesh; threads=threads)
    order_timeCS!(phi.v, mesh; threads=threads)
    order_timeCS!(phi.w, mesh; threads=threads)
    order_timeCS!(phi.p, mesh; threads=threads)

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            for j in 1:mesh.m1
                for k in 1:mesh.n1
                    if (i == mesh.l1)
                        @inbounds phi.fValues.uFaceTime3[i,j,k] = phi.fValues.uFaceTime2[i,j,k]
                        @inbounds phi.fValues.uFaceTime2[i,j,k] = phi.fValues.uFaceTime1[i,j,k]
                        @inbounds phi.fValues.uFaceTime1[i,j,k] = phi.fValues.uFace[i,j,k]
                        @inbounds phi.fValues.uFaceTime3[i+1,j,k] = phi.fValues.uFaceTime2[i+1,j,k]
                        @inbounds phi.fValues.uFaceTime2[i+1,j,k] = phi.fValues.uFaceTime1[i+1,j,k]
                        @inbounds phi.fValues.uFaceTime1[i+1,j,k] = phi.fValues.uFace[i+1,j,k]
                    else
                        @inbounds phi.fValues.uFaceTime3[i,j,k] = phi.fValues.uFaceTime2[i,j,k]
                        @inbounds phi.fValues.uFaceTime2[i,j,k] = phi.fValues.uFaceTime1[i,j,k]
                        @inbounds phi.fValues.uFaceTime1[i,j,k] = phi.fValues.uFace[i,j,k]
                    end
                    if (j == mesh.m1)
                        @inbounds phi.fValues.vFaceTime3[i,j,k] = phi.fValues.vFaceTime2[i,j,k]
                        @inbounds phi.fValues.vFaceTime2[i,j,k] = phi.fValues.vFaceTime1[i,j,k]
                        @inbounds phi.fValues.vFaceTime1[i,j,k] = phi.fValues.vFace[i,j,k]
                        @inbounds phi.fValues.vFaceTime3[i,j+1,k] = phi.fValues.vFaceTime2[i,j+1,k]
                        @inbounds phi.fValues.vFaceTime2[i,j+1,k] = phi.fValues.vFaceTime1[i,j+1,k]
                        @inbounds phi.fValues.vFaceTime1[i,j+1,k] = phi.fValues.vFace[i,j+1,k]
                    else
                        @inbounds phi.fValues.vFaceTime3[i,j,k] = phi.fValues.vFaceTime2[i,j,k]
                        @inbounds phi.fValues.vFaceTime2[i,j,k] = phi.fValues.vFaceTime1[i,j,k]
                        @inbounds phi.fValues.vFaceTime1[i,j,k] = phi.fValues.vFace[i,j,k]
                    end
                    if (k == mesh.n1)
                        @inbounds phi.fValues.wFaceTime3[i,j,k] = phi.fValues.wFaceTime2[i,j,k]
                        @inbounds phi.fValues.wFaceTime2[i,j,k] = phi.fValues.wFaceTime1[i,j,k]
                        @inbounds phi.fValues.wFaceTime1[i,j,k] = phi.fValues.wFace[i,j,k]
                        @inbounds phi.fValues.wFaceTime3[i,j,k+1] = phi.fValues.wFaceTime2[i,j,k+1]
                        @inbounds phi.fValues.wFaceTime2[i,j,k+1] = phi.fValues.wFaceTime1[i,j,k+1]
                        @inbounds phi.fValues.wFaceTime1[i,j,k+1] = phi.fValues.wFace[i,j,k+1]
                    else
                        @inbounds phi.fValues.wFaceTime3[i,j,k] = phi.fValues.wFaceTime2[i,j,k]
                        @inbounds phi.fValues.wFaceTime2[i,j,k] = phi.fValues.wFaceTime1[i,j,k]
                        @inbounds phi.fValues.wFaceTime1[i,j,k] = phi.fValues.wFace[i,j,k]
                    end
                end
            end
        end
    elseif !threads
        for i in 1:mesh.l1
            for j in 1:mesh.m1
                for k in 1:mesh.n1
                    if (i == mesh.l1)
                        @inbounds phi.fValues.uFaceTime3[i,j,k] = phi.fValues.uFaceTime2[i,j,k]
                        @inbounds phi.fValues.uFaceTime2[i,j,k] = phi.fValues.uFaceTime1[i,j,k]
                        @inbounds phi.fValues.uFaceTime1[i,j,k] = phi.fValues.uFace[i,j,k]
                        @inbounds phi.fValues.uFaceTime3[i+1,j,k] = phi.fValues.uFaceTime2[i+1,j,k]
                        @inbounds phi.fValues.uFaceTime2[i+1,j,k] = phi.fValues.uFaceTime1[i+1,j,k]
                        @inbounds phi.fValues.uFaceTime1[i+1,j,k] = phi.fValues.uFace[i+1,j,k]
                    else
                        @inbounds phi.fValues.uFaceTime3[i,j,k] = phi.fValues.uFaceTime2[i,j,k]
                        @inbounds phi.fValues.uFaceTime2[i,j,k] = phi.fValues.uFaceTime1[i,j,k]
                        @inbounds phi.fValues.uFaceTime1[i,j,k] = phi.fValues.uFace[i,j,k]
                    end
                    if (j == mesh.m1)
                        @inbounds phi.fValues.vFaceTime3[i,j,k] = phi.fValues.vFaceTime2[i,j,k]
                        @inbounds phi.fValues.vFaceTime2[i,j,k] = phi.fValues.vFaceTime1[i,j,k]
                        @inbounds phi.fValues.vFaceTime1[i,j,k] = phi.fValues.vFace[i,j,k]
                        @inbounds phi.fValues.vFaceTime3[i,j+1,k] = phi.fValues.vFaceTime2[i,j+1,k]
                        @inbounds phi.fValues.vFaceTime2[i,j+1,k] = phi.fValues.vFaceTime1[i,j+1,k]
                        @inbounds phi.fValues.vFaceTime1[i,j+1,k] = phi.fValues.vFace[i,j+1,k]
                    else
                        @inbounds phi.fValues.vFaceTime3[i,j,k] = phi.fValues.vFaceTime2[i,j,k]
                        @inbounds phi.fValues.vFaceTime2[i,j,k] = phi.fValues.vFaceTime1[i,j,k]
                        @inbounds phi.fValues.vFaceTime1[i,j,k] = phi.fValues.vFace[i,j,k]
                    end
                    if (k == mesh.n1)
                        @inbounds phi.fValues.wFaceTime3[i,j,k] = phi.fValues.wFaceTime2[i,j,k]
                        @inbounds phi.fValues.wFaceTime2[i,j,k] = phi.fValues.wFaceTime1[i,j,k]
                        @inbounds phi.fValues.wFaceTime1[i,j,k] = phi.fValues.wFace[i,j,k]
                        @inbounds phi.fValues.wFaceTime3[i,j,k+1] = phi.fValues.wFaceTime2[i,j,k+1]
                        @inbounds phi.fValues.wFaceTime2[i,j,k+1] = phi.fValues.wFaceTime1[i,j,k+1]
                        @inbounds phi.fValues.wFaceTime1[i,j,k+1] = phi.fValues.wFace[i,j,k+1]
                    else
                        @inbounds phi.fValues.wFaceTime3[i,j,k] = phi.fValues.wFaceTime2[i,j,k]
                        @inbounds phi.fValues.wFaceTime2[i,j,k] = phi.fValues.wFaceTime1[i,j,k]
                        @inbounds phi.fValues.wFaceTime1[i,j,k] = phi.fValues.wFace[i,j,k]
                    end
                end
            end
        end
    end

    return nothing
end
