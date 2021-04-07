"""

"""
function convergence_time(
    phi::UnionCSPhi,
    mesh::UnionCSMesh;
)
    maxvalue = maximum(abs.(phi.eval - phi.time1))

    return maxvalue
end
"""

"""
function convergence_time! end

function convergence_time!(
    phi::UnionCSPhi,
)
    phi.errorT = maximum(abs.(phi.eval - phi.time1))

    return nothing
end

function convergence_time!(
    velocity::CSVelocity1D;
)
    convergence_time!(velocity.u)

    return nothing
end

function convergence_time!(
    velocity::CSVelocity2D;
)
    convergence_time!(velocity.u)
    convergence_time!(velocity.v)

    return nothing
end

function convergence_time!(
    velocity::CSVelocity3D;
)
    convergence_time!(velocity.u)
    convergence_time!(velocity.v)
    convergence_time!(velocity.w)

    return nothing
end

function convergence_time!(
    args...;
)
    nt = length(args)

    for nn in 1:nt
        convergence_time!(args[nn])
    end

    return nothing
end
