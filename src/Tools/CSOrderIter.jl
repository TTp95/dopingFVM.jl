"""

"""
function order_iterCS! end

function order_iterCS!(
    phi::CSPhi1D,
    mesh::UnionCSMesh1D;
    threads = false,
)

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            @inbounds phi.iter[i] = phi.eval[i]
        end
    elseif !threads
        for i in 1:mesh.l1
            @inbounds phi.iter[i] = phi.eval[i]
        end
    end

    return nothing
end

function order_iterCS!(
    phi::CSPhi2D,
    mesh::UnionCSMesh2D;
    threads = false,
)

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            for j in 1:mesh.m1
                @inbounds phi.iter[i,j] = phi.eval[i,j]
            end
        end
    elseif !threads
        for i in 1:mesh.l1
            for j in 1:mesh.m1
                @inbounds phi.iter[i,j] = phi.eval[i,j]
            end
        end
    end

    return nothing
end

function order_iterCS!(
    phi::CSPhi3D,
    mesh::UnionCSMesh3D;
    threads = false,
)

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            for j in 1:mesh.m1
                for k in 1:mesh.n1
                    @inbounds phi.iter[i,j,k] = phi.eval[i,j,k]
                end
            end
        end
    elseif !threads
        for i in 1:mesh.l1
            for j in 1:mesh.m1
                for k in 1:mesh.n1
                    @inbounds phi.iter[i,j,k] = phi.eval[i,j,k]
                end
            end
        end
    end

    return nothing
end

function order_iterCS!(
    phi::CSVelocity1D,
    mesh::UnionCSMesh1D;
    threads = false,
)
    order_iterCS!(phi.u, mesh; threads=threads)
    order_iterCS!(phi.p, mesh; threads=threads)

    if threads
        Base.Threads.@threads for i in 1:(mesh.l1 + 1)
            @inbounds phi.fValues.uFaceIter[i] = phi.fValues.uFace[i]
        end
    elseif !threads
        for i in 1:(mesh.l1 + 1)
            @inbounds phi.fValues.uFaceIter[i] = phi.fValues.uFace[i]
        end
    end

    return nothing
end

function order_iterCS!(
    phi::CSVelocity2D,
    mesh::UnionCSMesh2D;
    threads = false,
)
    order_iterCS!(phi.u, mesh; threads=threads)
    order_iterCS!(phi.v, mesh; threads=threads)
    order_iterCS!(phi.p, mesh; threads=threads)

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            for j in 1:mesh.m1
                if (i == mesh.l1)
                    @inbounds phi.fValues.uFaceIter[i,j] = phi.fValues.uFace[i,j]
                    @inbounds phi.fValues.uFaceIter[i+1,j] = phi.fValues.uFace[i+1,j]
                else
                    @inbounds phi.fValues.uFaceIter[i,j] = phi.fValues.uFace[i,j]
                end
                if (j == mesh.m1)
                    @inbounds phi.fValues.vFaceIter[i,j] = phi.fValues.vFace[i,j]
                    @inbounds phi.fValues.vFaceIter[i,j+1] = phi.fValues.vFace[i,j+1]
                else
                    @inbounds phi.fValues.vFaceIter[i,j] = phi.fValues.vFace[i,j]
                end
            end
        end
    elseif !threads
        for i in 1:mesh.l1
            for j in 1:mesh.m1
                if (i == mesh.l1)
                    @inbounds phi.fValues.uFaceIter[i,j] = phi.fValues.uFace[i,j]
                    @inbounds phi.fValues.uFaceIter[i+1,j] = phi.fValues.uFace[i+1,j]
                else
                    @inbounds phi.fValues.uFaceIter[i,j] = phi.fValues.uFace[i,j]
                end
                if (j == mesh.m1)
                    @inbounds phi.fValues.vFaceIter[i,j] = phi.fValues.vFace[i,j]
                    @inbounds phi.fValues.vFaceIter[i,j+1] = phi.fValues.vFace[i,j+1]
                else
                    @inbounds phi.fValues.vFaceIter[i,j] = phi.fValues.vFace[i,j]
                end
            end
        end
    end

    return nothing
end

function order_iterCS!(
    phi::CSVelocity3D,
    mesh::UnionCSMesh3D;
    threads = false,
)
    order_iterCS!(phi.u, mesh; threads=threads)
    order_iterCS!(phi.v, mesh; threads=threads)
    order_iterCS!(phi.w, mesh; threads=threads)
    order_iterCS!(phi.p, mesh; threads=threads)

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            for j in 1:mesh.m1
                for k in 1:mesh.n1
                    if (i == mesh.l1)
                        @inbounds phi.fValues.uFaceIter[i,j,k] = phi.fValues.uFace[i,j,k]
                        @inbounds phi.fValues.uFaceIter[i+1,j,k] = phi.fValues.uFace[i+1,j,k]
                    else
                        @inbounds phi.fValues.uFaceIter[i,j,k] = phi.fValues.uFace[i,j,k]
                    end
                    if (j == mesh.m1)
                        @inbounds phi.fValues.vFaceIter[i,j,k] = phi.fValues.vFace[i,j,k]
                        @inbounds phi.fValues.vFaceIter[i,j+1,k] = phi.fValues.vFace[i,j+1,k]
                    else
                        @inbounds phi.fValues.vFaceIter[i,j,k] = phi.fValues.vFace[i,j,k]
                    end
                    if (k == mesh.n1)
                        @inbounds phi.fValues.wFaceIter[i,j,k] = phi.fValues.wFace[i,j,k]
                        @inbounds phi.fValues.wFaceIter[i,j,k+1] = phi.fValues.wFace[i,j,k+1]
                    else
                        @inbounds phi.fValues.wFaceIter[i,j,k] = phi.fValues.wFace[i,j,k]
                    end
                end
            end
        end
    elseif !threads
        for i in 1:mesh.l1
            for j in 1:mesh.m1
                for k in 1:mesh.n1
                    if (i == mesh.l1)
                        @inbounds phi.fValues.uFaceIter[i,j,k] = phi.fValues.uFace[i,j,k]
                        @inbounds phi.fValues.uFaceIter[i+1,j,k] = phi.fValues.uFace[i+1,j,k]
                    else
                        @inbounds phi.fValues.uFaceIter[i,j,k] = phi.fValues.uFace[i,j,k]
                    end
                    if (j == mesh.m1)
                        @inbounds phi.fValues.vFaceIter[i,j,k] = phi.fValues.vFace[i,j,k]
                        @inbounds phi.fValues.vFaceIter[i,j+1,k] = phi.fValues.vFace[i,j+1,k]
                    else
                        @inbounds phi.fValues.vFaceIter[i,j,k] = phi.fValues.vFace[i,j,k]
                    end
                    if (k == mesh.n1)
                        @inbounds phi.fValues.wFaceIter[i,j,k] = phi.fValues.wFace[i,j,k]
                        @inbounds phi.fValues.wFaceIter[i,j,k+1] = phi.fValues.wFace[i,j,k+1]
                    else
                        @inbounds phi.fValues.wFaceIter[i,j,k] = phi.fValues.wFace[i,j,k]
                    end
                end
            end
        end
    end

    return nothing
end
