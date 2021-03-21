"""

"""
function convergence_time(
    phi::UnionCSPhi,
    mesh::UnionCSMesh;
)
    maxvalue = maximum(abs.(phi.eval - phi.time1))

    return maxvalue
end
