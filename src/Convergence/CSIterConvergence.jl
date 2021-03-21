"""
convergence_iter(phi::Union{CSPhi1D,CSPhi2D,CSPhi3D}, mesh::Union{CSMesh1D,CSMesh1DImmutable})

"""
function convergence_iter(
    phi::UnionCSPhi,
    mesh::UnionCSMesh;
)
    maxValue = maximum(abs.(phi.eval - phi.iter))

    return maxValue
end

"""
convergence_relative_iter(phi::Union{CSPhi1D,CSPhi2D,CSPhi3D}, mesh::Union{CSMesh1D,CSMesh1DImmutable})

"""
function convergence_relative_iter(
    phi::UnionCSPhi,
    mesh::UnionCSMesh;
)
    maxValue = maximum(abs.((phi.eval - phi.iter)./(phi.eval + 1.0e-10)))

    return maxValue
end
