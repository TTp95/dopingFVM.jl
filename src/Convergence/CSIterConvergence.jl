"""

"""
function convergence_iter(
    phi::UnionCSPhi,
)
    maxValue = maximum(abs.(phi.eval - phi.iter))

    return maxValue
end

"""

"""
function convergence_iter! end

function convergence_iter!(
    phi::UnionCSPhi,
)
    phi.errorI = maximum(abs.(phi.eval - phi.iter))

    return nothing
end

function convergence_iter!(
    velocity::CSVelocity1D;
)
    convergence_iter!(velocity.u)

    return nothing
end

function convergence_iter!(
    velocity::CSVelocity2D;
)
    convergence_iter!(velocity.u)
    convergence_iter!(velocity.v)

    return nothing
end

function convergence_iter!(
    velocity::CSVelocity3D;
)
    convergence_iter!(velocity.u)
    convergence_iter!(velocity.v)
    convergence_iter!(velocity.w)

    return nothing
end

function convergence_iter!(
    args...;
)
    nt = length(args)

    for nn in 1:nt
        convergence_iter!(args[nn])
    end

    return nothing
end

"""

"""
function convergence_relativeIter(
    phi::UnionCSPhi,
)
    maxValue = maximum(abs.((phi.eval - phi.iter)./(phi.eval + 1.0e-10)))

    return maxValue
end

"""

"""
function convergence_relativeIter! end

function convergence_relativeIter!(
    phi::UnionCSPhi,
)
    phi.errorR = maximum(abs.((phi.eval - phi.iter)./(phi.eval + 1.0e-10)))

    return nothing
end

function convergence_relativeIter!(
    velocity::CSVelocity1D;
)
    convergence_relativeIter!(velocity.u)

    return nothing
end

function convergence_relativeIter!(
    velocity::CSVelocity2D;
)
    convergence_relativeIter!(velocity.u)
    convergence_relativeIter!(velocity.v)

    return nothing
end

function convergence_relativeIter!(
    velocity::CSVelocity3D;
)
    convergence_relativeIter!(velocity.u)
    convergence_relativeIter!(velocity.v)
    convergence_relativeIter!(velocity.w)

    return nothing
end

function convergence_relativeIter!(
    args...;
)
    nt = length(args)

    for nn in 1:nt
        convergence_relativeIter!(args[nn])
    end

    return nothing
end
