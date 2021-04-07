"""

"""
function check_iterConvergence end

function check_iterConvergence(
    phi::UnionCSPhi;
    maxvalue::AbstractFloat = phi.errorI,
)
    if (maxvalue <= phi.convergence)
        solve = true
    else
        solve = false
    end

    return solve
end

function check_iterConvergence(
    velocity::CSVelocity1D;
    mass::AbstractFloat = -1.0,
)
    check1 = check_iterConvergence(velocity.u)
    check0 = true

    if (mass > 0.0)
        check0 = (mass <= velocity.p.convergence) ? true : false
    end

    if check1 && check0
        solve = true
    else
        solve = false
    end

    return solve
end

function check_iterConvergence(
    velocity::CSVelocity2D;
    mass::AbstractFloat = -1.0,
)
    check1 = check_iterConvergence(velocity.u)
    check2 = check_iterConvergence(velocity.v)
    check0 = true

    if (mass > 0.0)
        check0 = (mass <= velocity.p.convergence) ? true : false
    end

    if check0 && check1 && check2
        solve = true
    else
        solve = false
    end

    return solve
end

function check_iterConvergence(
    velocity::CSVelocity3D;
    mass::AbstractFloat = -1.0,
)
    check1 = check_iterConvergence(velocity.u)
    check2 = check_iterConvergence(velocity.v)
    check3 = check_iterConvergence(velocity.w)
    check0 = true

    if (mass > 0.0)
        check0 = (mass <= velocity.p.convergence) ? true : false
    end

    if check0 && check1 && check2 && check3
        solve = true
    else
        solve = false
    end

    return solve
end

"""

"""
function check_iterConvergence!(
    system::SystemControl,
    args...;
    mass::AbstractFloat = -1.0,
)
    nt = length(args)

    solve = true

    for nn in 1:nt
        if (mass > 0.0) && (typeof(args[nn]) <: UnionCSVelocity)
            bool_aux = check_iterConvergence(args[nn]; mass = mass)
            solve *= bool_aux
        else
            bool_aux = check_iterConvergence(args[nn])
            solve *= bool_aux
        end
    end

    system.solvedTimeStep = solve

    return nothing
end

"""

"""
function check_timeConvergence end

function check_timeConvergence(
    phi::UnionCSPhi;
    maxvalue::AbstractFloat = phi.errorT,
)
    if (maxvalue <= phi.convergenceTime)
        solve = true
    else
        solve = false
    end

    return solve
end

function check_timeConvergence(
    velocity::CSVelocity1D;
    mass::AbstractFloat = -1.0,
)
    solve = check_timeConvergence(velocity.u)

    return solve
end

function check_timeConvergence(
    velocity::CSVelocity2D;
    mass::AbstractFloat = -1.0,
)
    check1 = check_timeConvergence(velocity.u)
    check2 = check_timeConvergence(velocity.v)

    if check1 && check2
        solve = true
    else
        solve = false
    end

    return solve
end

function check_timeConvergence(
    velocity::CSVelocity3D;
    mass::AbstractFloat = -1.0,
)
    check1 = check_timeConvergence(velocity.u)
    check2 = check_timeConvergence(velocity.v)
    check3 = check_timeConvergence(velocity.w)

    if check1 && check2 && check3
        solve = true
    else
        solve = false
    end

    return solve
end

"""

"""
function check_timeConvergence!(
    system::SystemControl,
    args...;
)
    nt = length(args)

    solve = true

    for nn in 1:nt
        bool_aux = check_timeConvergence(args[nn])
        solve *= bool_aux
    end

    system.solvedProblem = solve

    return nothing
end

"""

"""
function check_relativeConvergence end

function check_relativeConvergence(
    phi::UnionCSPhi;
    maxvalue::AbstractFloat = phi.errorR,
)
    if (maxvalue <= phi.convergenceRelative)
        solve = true
    else
        solve = false
    end

    return solve
end

function check_relativeConvergence(
    velocity::CSVelocity1D;
    mass::AbstractFloat = -1.0,
)
    check1 = check_relativeConvergence(velocity.u)
    check0 = true

    if (mass > 0.0)
        check0 = (mass <= velocity.p.convergence) ? true : false
    end

    if check1 && check0
        solve = true
    else
        solve = false
    end

    return solve
end

function check_relativeConvergence(
    velocity::CSVelocity2D;
    mass::AbstractFloat = -1.0,
)
    check1 = check_relativeConvergence(velocity.u)
    check2 = check_relativeConvergence(velocity.v)
    check0 = true

    if (mass > 0.0)
        check0 = (mass <= velocity.p.convergence) ? true : false
    end

    if check0 && check1 && check2
        solve = true
    else
        solve = false
    end

    return solve
end

function check_relativeConvergence(
    velocity::CSVelocity3D;
    mass::AbstractFloat = -1.0,
)
    check1 = check_relativeConvergence(velocity.u)
    check2 = check_relativeConvergence(velocity.v)
    check3 = check_relativeConvergence(velocity.w)
    check0 = true

    if (mass > 0.0)
        check0 = (mass <= velocity.p.convergence) ? true : false
    end

    if check0 && check1 && check2 && check3
        solve = true
    else
        solve = false
    end

    return solve
end

"""

"""
function check_relativeConvergence!(
    system::SystemControl,
    args...;
    mass::AbstractFloat = -1.0,
)
    nt = length(args)

    solve = true

    for nn in 1:nt
        if (mass > 0.0) && (typeof(args[nn]) <: UnionCSVelocity)
            bool_aux = check_relativeConvergence(args[nn]; mass = mass)
            solve *= bool_aux
        else
            bool_aux = check_relativeConvergence(args[nn])
            solve *= bool_aux
        end
    end

    system.solvedTimeStep = solve

    return nothing
end
