"""

"""
@inline function _convection_upwind_neighbors_(
    mflux::AbstractFloat,
)
    aF =  -1.0 * max((-1.0 * mflux), 0.0)
    aFc =  max(mflux, 0.0)
    b = 0.0

    return aF, aFc, b
end

"""

"""
function _convection_upwind_central_ end

function _convection_upwind_central_(
    i::Signed,
    velocityU::Array{<:AbstractFloat,1},
    phi::CSPhi1D,
    bounds::Dict{String,BoundsStructured},
    mesh::UnionCSMesh1D;
    T::Type{<:AbstractFloat} = Float64,
)
    nbounds = phi.nbounds[i]
    ax = zeros(T, nbounds)
    bx = zeros(T, nbounds)

    for nn in  1:nbounds
        @inbounds ax[nn], bx[nn] =  _convection_upwind_bounds_(
            i,
            velocityU,
            phi,
            bounds["$(phi.gIndex[i])g$(nn)"],
            mesh,
        )
    end

    aC = sum(ax)
    b = sum(bx)

    return aC, b
end


function _convection_upwind_central_(
    i::Signed,
    j::Signed,
    velocityU::Array{<:AbstractFloat,2},
    velocityV::Array{<:AbstractFloat,2},
    phi::CSPhi2D,
    bounds::Dict{String,BoundsStructured},
    mesh::UnionCSMesh2D;
    T::Type{<:AbstractFloat} = Float64,
)
    nbounds = phi.nbounds[i,j]
    ax = zeros(T, nbounds)
    bx = zeros(T, nbounds)

    for nn in  1:nbounds
        @inbounds ax[nn], bx[nn] =  _convection_upwind_bounds_(
            i,
            j,
            velocityU,
            velocityV,
            phi,
            bounds["$(phi.gIndex[i,j])g$(nn)"],
            mesh,
        )
    end

    aC = sum(ax)
    b = sum(bx)

    return aC, b
end

function _convection_upwind_central_(
    i::Signed,
    j::Signed,
    k::Signed,
    velocityU::Array{<:AbstractFloat,3},
    velocityV::Array{<:AbstractFloat,3},
    velocityW::Array{<:AbstractFloat,3},
    phi::CSPhi3D,
    bounds::Dict{String,BoundsStructured},
    mesh::UnionCSMesh3D;
    T::Type{<:AbstractFloat} = Float64,
)
    nbounds = phi.nbounds[i,j,k]
    ax = zeros(T, nbounds)
    bx = zeros(T, nbounds)

    for nn in  1:nbounds
        @inbounds ax[nn], bx[nn] =  _convection_upwind_bounds_(
            i,
            j,
            k,
            velocityU,
            velocityV,
            velocityW,
            phi,
            bounds["$(phi.gIndex[i,j,k])g$(nn)"],
            mesh,
        )
    end

    aC = sum(ax)
    b = sum(bx)

    return aC, b
end
