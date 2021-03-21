"""

"""
function _diffusion_centralDifference_neighbors_(
    gamma::AbstractFloat,
    lenghts::Union{Array{<:AbstractFloat,1}, Array{<:AbstractFloat,2}},
    area::AbstractFloat,
)
    dCF = 0.5 * (lenghts[1] + lenghts[2])
    aF = -1.0 * gamma * (area / dCF)
    b = 0.0

    return aF, b
end

"""

"""
function _diffusion_centralDifference_central_ end

function _diffusion_centralDifference_central_(
    i::Signed,
    phi::CSPhi1D,
    bounds::Dict{String,BoundsStructured},
    mesh::UnionCSMesh1D;
    T::Type{<:AbstractFloat} = Float64,
)
    nbounds = phi.nbounds[i]
    ax = zeros(T, nbounds)
    bx = zeros(T, nbounds)

    for nn in  1:nbounds
        @inbounds ax[nn], bx[nn] =  _diffusion_centralDifference_bounds_(
            i,
            phi,
            bounds["$(phi.gIndex[i])g$(nn)"],
            mesh,
        )
    end

    aC = sum(ax)
    b = sum(bx)

    return aC, b
end

function _diffusion_centralDifference_central_(
    i::Signed,
    j::Signed,
    phi::CSPhi2D,
    bounds::Dict{String,BoundsStructured},
    mesh::UnionCSMesh2D;
    T::Type{<:AbstractFloat} = Float64,
)
    nbounds = phi.nbounds[i,j]
    ax = zeros(T, nbounds)
    bx = zeros(T, nbounds)

    for nn in  1:nbounds
        @inbounds ax[nn], bx[nn] =  _diffusion_centralDifference_bounds_(
            i,
            j,
            phi,
            bounds["$(phi.gIndex[i,j])g$(nn)"],
            mesh,
        )
    end

    aC = sum(ax)
    b = sum(bx)

    return aC, b
end

function _diffusion_centralDifference_central_(
    i::Signed,
    j::Signed,
    k::Signed,
    phi::CSPhi3D,
    bounds::Dict{String,BoundsStructured},
    mesh::UnionCSMesh3D;
    T::Type{<:AbstractFloat} = Float64,
)
    nbounds = phi.nbounds[i,j,k]
    ax = zeros(T, nbounds)
    bx = zeros(T, nbounds)

    for nn in  1:nbounds
        @inbounds ax[nn], bx[nn] =  _diffusion_centralDifference_bounds_(
            i,
            j,
            k,
            phi,
            bounds["$(phi.gIndex[i,j,k])g$(nn)"],
            mesh,
        )
    end

    aC = sum(ax)
    b = sum(bx)

    return aC, b
end
