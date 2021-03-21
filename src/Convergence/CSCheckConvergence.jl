"""

"""
function check_iterConvergence(
    maxvalue::AbstractFloat,
    phi::UnionCSPhi,
    mesh::UnionCSMesh;
)
    if (maxvalue <= phi.convergence)
        solve = true
    else
        solve = false
    end

    return solve
end

"""

"""
function check_timeConvergence(
    maxvalue::AbstractFloat,
    phi::UnionCSPhi,
    mesh::UnionCSMesh;
)
    if (maxvalue <= phi.convergenceTime)
        solve = true
    else
        solve = false
    end

    return solve
end

"""

"""
function check_timeConvergenceRelavite(
    maxvalue::AbstractFloat,
    phi::UnionCSPhi,
    mesh::UnionCSMesh;
)
    if (maxvalue <= phi.convergenceRelative)
        solve = true
    else
        solve = false
    end

    return solve
end
