"""

"""
function check_bounds! end

function check_bounds!(
    phi::CSPhi1D,
    mesh::UnionCSMesh1D,
    dictBounds::Dict{String,BoundsStructured};
    threads = false,
)

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            if phi.bounds[i]
                for nn in 1:phi.nbounds[i]
                    id = phi.gIndex[i]
                    _evaluate_bounds!_(
                        i, phi, mesh, dictBounds["$(id)g$(nn)"];
                    )
                end
            end
        end
    elseif !threads
        for i in 1:mesh.l1
            if phi.bounds[i]
                for nn in 1:phi.nbounds[i]
                    id = phi.gIndex[i]
                    _evaluate_bounds!_(
                        i, phi, mesh, dictBounds["$(id)g$(nn)"];
                    )
                end
            end
        end
    end

    return nothing
end

function check_bounds!(
    phi::CSPhi2D,
    mesh::UnionCSMesh2D,
    dictBounds::Dict{String,BoundsStructured};
    threads = false,
)

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            for j in 1:mesh.m1
                if phi.bounds[i,j]
                    for nn in 1:phi.nbounds[i,j]
                        id = phi.gIndex[i,j]
                        _evaluate_bounds!_(
                            i, j, phi, mesh, dictBounds["$(id)g$(nn)"];
                        )
                    end
                end
            end
        end
    elseif !threads
        for i in 1:mesh.l1
            for j in 1:mesh.m1
                if phi.bounds[i,j]
                    for nn in 1:phi.nbounds[i,j]
                        id = phi.gIndex[i,j]
                        _evaluate_bounds!_(
                            i, j, phi, mesh, dictBounds["$(id)g$(nn)"];
                        )
                    end
                end
            end
        end
    end

    return nothing
end

function check_bounds!(
    phi::CSPhi3D,
    mesh::UnionCSMesh3D,
    dictBounds::Dict{String,BoundsStructured};
    threads = false,
)

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            for j in 1:mesh.m1
                for j in 1:mesh.n1
                    if phi.bounds[i,j,k]
                        for nn in 1:phi.nbounds[i,j,k]
                            id = phi.gIndex[i,j,k]
                            _evaluate_bounds!_(
                                i, j, k, phi, mesh, dictBounds["$(id)g$(nn)"];
                            )
                        end
                    end
                end
            end
        end
    elseif !threads
        for i in 1:mesh.l1
            for j in 1:mesh.m1
                for j in 1:mesh.n1
                    if phi.bounds[i,j,k]
                        for nn in 1:phi.nbounds[i,j,k]
                            id = phi.gIndex[i,j,k]
                            _evaluate_bounds!_(
                                i, j, k, phi, mesh, dictBounds["$(id)g$(nn)"];
                            )
                        end
                    end
                end
            end
        end
    end

    return nothing
end
