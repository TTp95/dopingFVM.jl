"""

"""
function assign_convergenceParameters!(
    phi::UnionCSPhi,
    args...,
)
    nt = length(args)

    if (nt > 4)
        nt = 4
    end

    if (nt == 1)
        phi.convergence = args[1]
    elseif (nt == 2)
        phi.convergence = args[1]
        phi.convergenceRelative = args[2]
    elseif (nt == 3)
        phi.convergence = args[1]
        phi.convergenceRelative = args[2]
        phi.convergenceTime = args[3]
    elseif (nt == 3)
        phi.convergence = args[1]
        phi.convergenceRelative = args[2]
        phi.convergenceTime = args[3]
        phi.tolerance = args[4]
    end

    return nothing
end
