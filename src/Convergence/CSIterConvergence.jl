"""
convergence_iterCS(phi::Union{CSPhi1D,CSPhi2D,CSPhi3D}, mesh::Union{CSMesh1D,CSMesh1DImmutable})

"""
function convergence_iterCS(
    phi::UnionCSPhi,
    mesh::UnionCSMesh;
)
    maxValue = maximum(abs.(phi.eval - phi.iter))

    return maxValue
end

"""
convergence_relative_iterCS(phi::Union{CSPhi1D,CSPhi2D,CSPhi3D}, mesh::Union{CSMesh1D,CSMesh1DImmutable})

"""
function convergence_relative_iterCS(
    phi::UnionCSPhi,
    mesh::UnionCSMesh;
)
    maxValue = maximum(abs.((phi.eval - phi.iter)./(phi.eval + 1.0e-10)))

    return maxValue
end
