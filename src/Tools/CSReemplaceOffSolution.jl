"""

"""
function remplace_offSolution! end

function remplace_offSolution!(
    phi::CSPhi1D,
    mesh::UnionCSMesh1D,
    offNodeValue::AbstractFloat;
    threads = false,
)

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            if !phi.onoff
                phi.eval[i] = offNodeValue
            end
        end
    elseif !threads
        for i in 1:mesh.l1
            if !phi.onoff
                phi.eval[i] = offNodeValue
            end
        end
    end

    return nothing
end

function remplace_offSolution!(
    phi::CSPhi2D,
    mesh::UnionCSMesh2D,
    offNodeValue::AbstractFloat;
    threads = false,
)

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            for j in 1:mesh.m1
                if !phi.onoff
                    phi.eval[i,j] = offNodeValue
                end
            end
        end
    elseif !threads
        for i in 1:mesh.l1
            for j in 1:mesh.m1
                if !phi.onoff
                    phi.eval[i,j] = offNodeValue
                end
            end
        end
    end

    return nothing
end

function remplace_offSolution!(
    phi::CSPhi3D,
    mesh::UnionCSMesh3D,
    offNodeValue::AbstractFloat;
    threads = false,
)

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            for j in 1:mesh.m1
                for k in 1:mesh.n1
                    if !phi.onoff
                        phi.eval[i,j,k] = offNodeValue
                    end
                end
            end
        end
    elseif !threads
        for i in 1:mesh.l1
            for j in 1:mesh.m1
                for k in 1:mesh.n1
                    if !phi.onoff
                        phi.eval[i,j,k] = offNodeValue
                    end
                end
            end
        end
    end

    return nothing
end

function remplace_offSolution!(
    phi::CSPhi1D,
    mesh::UnionCSMesh1D,
    offNodeValue::Array{AbstractFloat,1};
    threads = false,
)

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            if !phi.onoff
                phi.eval[i] = offNodeValue[i]
            end
        end
    elseif !threads
        for i in 1:mesh.l1
            if !phi.onoff
                phi.eval[i] = offNodeValue[i]
            end
        end
    end

    return nothing
end

function remplace_offSolution!(
    phi::CSPhi2D,
    mesh::UnionCSMesh2D,
    offNodeValue::Array{AbstractFloat,2};
    threads = false,
)

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            for j in 1:mesh.m1
                if !phi.onoff
                    phi.eval[i,j] = offNodeValue[i,j]
                end
            end
        end
    elseif !threads
        for i in 1:mesh.l1
            for j in 1:mesh.m1
                if !phi.onoff
                    phi.eval[i,j] = offNodeValue[i,j]
                end
            end
        end
    end

    return nothing
end

function remplace_offSolution!(
    phi::CSPhi3D,
    mesh::UnionCSMesh3D,
    offNodeValue::Array{AbstractFloat,3};
    threads = false,
)

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            for j in 1:mesh.m1
                for k in 1:mesh.n1
                    if !phi.onoff
                        phi.eval[i,j,k] = offNodeValue[i,j,k]
                    end
                end
            end
        end
    elseif !threads
        for i in 1:mesh.l1
            for j in 1:mesh.m1
                for k in 1:mesh.n1
                    if !phi.onoff
                        phi.eval[i,j,k] = offNodeValue[i,j,k]
                    end
                end
            end
        end
    end

    return nothing
end
