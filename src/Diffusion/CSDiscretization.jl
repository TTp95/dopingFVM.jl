"""

"""
function discretize_diffusion end

function discretize_diffusion(
    phi::UnionCSPhi,
    bounds::Dict{String,BoundsStructured},
    material::UnionCSMaterialAll,
    mesh::UnionCSMesh;
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
    mthreads::Bool = false,
    sparrays::Bool = true,
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
            mthreads = mthreads,
            sparrays = sparrays,
            interpolation = interpolation,
        )
    else
        error("Diffusion scheme number $(scheme) unimplemented.")
    end

    return A, b
end
