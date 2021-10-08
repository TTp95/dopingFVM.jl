"""

"""
function discretize_PPC_pressureEquation end

function discretize_PPC_pressureEquation(
    velocity::CSVelocity1D,
    mesh::UnionCSMesh1D,
    deltat::DeltaTime,
    material::UnionCSMaterial;
    velocityU::Array{<:AbstractFloat,1} = velocity.fValues.uFace,
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
    transientScheme::Signed = 1,
    interpolation::Signed = 1,
)
    A = _projection_PC_pressure_laplacian_(
        velocity.p,
        mesh;
        T = T,
        N = N,
    )

    b = _projection_PC_velocity_divergence_(
        velocity,
        mesh,
        deltat,
        material;
        velocityU = velocityU,
        T = T,
        transientScheme = transientScheme,
        interpolation = interpolation,
    )

    return A, b
end

function discretize_PPC_pressureEquation(
    velocity::CSVelocity2D,
    mesh::UnionCSMesh2D,
    deltat::DeltaTime,
    material::UnionCSMaterial;
    velocityU::Array{<:AbstractFloat,2} = velocity.fValues.uFace,
    velocityV::Array{<:AbstractFloat,2} = velocity.fValues.vFace,
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
    transientScheme::Signed = 1,
    interpolation::Signed = 1,
)
    A = _projection_PC_pressure_laplacian_(
        velocity.p,
        mesh;
        T = T,
        N = N,
    )

    b = _projection_PC_velocity_divergence_(
        velocity,
        mesh,
        deltat,
        material;
        velocityU = velocityU,
        velocityV = velocityV,
        T = T,
        transientScheme = transientScheme,
        interpolation = interpolation,
    )

    return A, b
end

function discretize_PPC_pressureEquation(
    velocity::CSVelocity3D,
    mesh::UnionCSMesh3D,
    deltat::DeltaTime,
    material::UnionCSMaterial;
    velocityU::Array{<:AbstractFloat,3} = velocity.fValues.uFace,
    velocityV::Array{<:AbstractFloat,3} = velocity.fValues.vFace,
    velocityW::Array{<:AbstractFloat,3} = velocity.fValues.wFace,
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
    transientScheme::Signed = 1,
    interpolation::Signed = 1,
)
    A = _projection_PC_pressure_laplacian_(
        velocity.p,
        mesh;
        T = T,
        N = N,
    )

    b = _projection_PC_velocity_divergence_(
        velocity,
        mesh,
        deltat,
        material;
        velocityU = velocityU,
        velocityV = velocityV,
        velocityW = velocityW,
        T = T,
        transientScheme = transientScheme,
        interpolation = interpolation,
    )

    return A, b
end
