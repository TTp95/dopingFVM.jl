"""

"""
function discretize_diffusion end

function discretize_diffusion(
    phi::UnionCSPhi,
    bounds::Dict{String,BoundsStructured},
    material::UnionCSMaterial,
    mesh::UnionCSMesh;
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
    scheme::Signed = 1,
    interpolation::Signed = 2,
)
    if (scheme == 1)
        A, b = _discretize_diffusion_centralDifference_(
            phi,
            bounds,
            material,
            mesh,;
            T = T,
            N = N,
            interpolation = interpolation,
        )
    else
        error("Diffusion scheme number $(scheme) unimplemented.")
    end

    return A, b
end
