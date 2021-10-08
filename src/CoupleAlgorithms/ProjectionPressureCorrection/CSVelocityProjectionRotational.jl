"""

"""
function velocityProjection_PPC_Rotational! end

function velocityProjection_PPC_Rotational!(
    field::Array{<:AbstractFloat,1},
    velocity::CSVelocity1D,
    mesh::UnionCSMesh1D,
    deltat::DeltaTime,
    material::UnionCSMaterial;
    pressure::AbstractArray = velocity.p.time1,
    velocityU::Array{<:AbstractFloat,1} = velocity.fValues.uFace,
    T::Type{<:AbstractFloat} = Float64,
    transientScheme::Signed = 1,
    interpolation::Signed = 2,
)
    if (transientScheme == 1)
        coef = deltat.dt1
    elseif (transientScheme == 2)
        coef = ((deltat.dt1 * (deltat.dt1 + deltat.dt2)) / deltat.dt2)
    elseif (transientScheme == 3)
        coef = 1.0 / ((1.0 / deltat.dt1 ) + ( 1.0 / (deltat.dt1 + deltat.dt2)))
    else
        error("Transient scheme unimplemented...")
    end

    array_field = zeros(T, mesh.l1)
    array_gradFieldx = zeros(T, mesh.l1)

    vector_to_phi!(field, velocity.p, mesh; phisolution = array_field)

    vector_fieldx = pressure_phi_gradient(velocity.p, mesh; phisolution = array_field, T = T)

    vector_to_phi!(vector_fieldx, velocity.p, mesh; phisolution = array_gradFieldx)

    array_div_vel = divergence_velocityToArray(
        velocity,
        mesh,
        material;
        velocityU = velocityU,
        T = T,
        interpolation = interpolation,
    )

    # Pressure
    velocity.p.eval .= array_field + pressure - array_div_vel

    # Velocity Projection
    velocity.u.eval .= velocity.u.eval - (coef * ((1.0 ./ material.ρ) .* array_gradFieldx) ./ mesh.vol)

    return nothing
end

function velocityProjection_PPC_Rotational!(
    field::Array{<:AbstractFloat,1},
    velocity::CSVelocity2D,
    mesh::UnionCSMesh2D,
    deltat::DeltaTime,
    material::UnionCSMaterial;
    pressure::AbstractArray = velocity.p.time1,
    velocityU::Array{<:AbstractFloat,2} = velocity.fValues.uFace,
    velocityV::Array{<:AbstractFloat,2} = velocity.fValues.vFace,
    T::Type{<:AbstractFloat} = Float64,
    transientScheme::Signed = 1,
    interpolation::Signed = 2,
)
    if (transientScheme == 1)
        coef = deltat.dt1
    elseif (transientScheme == 2)
        coef = ((deltat.dt1 * (deltat.dt1 + deltat.dt2)) / deltat.dt2)
    elseif (transientScheme == 3)
        coef = 1.0 / ((1.0 / deltat.dt1 ) + ( 1.0 / (deltat.dt1 + deltat.dt2)))
    else
        error("Transient scheme unimplemented...")
    end

    array_field = zeros(T, mesh.l1, mesh.m1)
    array_gradFieldx = zeros(T, mesh.l1, mesh.m1)
    array_gradFieldy = zeros(T, mesh.l1, mesh.m1)

    vector_to_phi!(field, velocity.p, mesh; phisolution = array_field)

    vector_fieldx, vector_fieldy = pressure_phi_gradient(velocity.p, mesh; phisolution = array_field, T = T)

    vector_to_phi!(vector_fieldx, velocity.p, mesh; phisolution = array_gradFieldx)
    vector_to_phi!(vector_fieldy, velocity.p, mesh; phisolution = array_gradFieldy)

    array_div_vel = divergence_velocityToArray(
        velocity,
        mesh,
        material;
        velocityU = velocityU,
        velocityV = velocityV,
        T = T,
        interpolation = interpolation,
    )

    # Pressure
    velocity.p.eval .= array_field + pressure - array_div_vel

    # Velocity Projection
    velocity.u.eval .= velocity.u.eval - (coef * ((1.0 ./ material.ρ) .* array_gradFieldx) ./ mesh.vol)
    velocity.v.eval .= velocity.v.eval - (coef * ((1.0 ./ material.ρ) .* array_gradFieldy) ./ mesh.vol)

    return nothing
end

function velocityProjection_PPC_Rotational!(
    field::Array{<:AbstractFloat,1},
    velocity::CSVelocity3D,
    mesh::UnionCSMesh3D,
    deltat::DeltaTime,
    material::UnionCSMaterial;
    pressure::AbstractArray = velocity.p.time1,
    velocityU::Array{<:AbstractFloat,3} = velocity.fValues.uFace,
    velocityV::Array{<:AbstractFloat,3} = velocity.fValues.vFace,
    velocityW::Array{<:AbstractFloat,3} = velocity.fValues.wFace,
    T::Type{<:AbstractFloat} = Float64,
    transientScheme::Signed = 1,
    interpolation::Signed = 2,
)
    if (transientScheme == 1)
        coef = deltat.dt1
    elseif (transientScheme == 2)
        coef = ((deltat.dt1 * (deltat.dt1 + deltat.dt2)) / deltat.dt2)
    elseif (transientScheme == 3)
        coef = 1.0 / ((1.0 / deltat.dt1 ) + ( 1.0 / (deltat.dt1 + deltat.dt2)))
    else
        error("Transient scheme unimplemented...")
    end

    array_field = zeros(T, mesh.l1, mesh.m1, mesh.n1)
    array_gradFieldx = zeros(T, mesh.l1, mesh.m1, mesh.n1)
    array_gradFieldy = zeros(T, mesh.l1, mesh.m1, mesh.n1)
    array_gradFieldz = zeros(T, mesh.l1, mesh.m1, mesh.n1)

    vector_to_phi!(field, velocity.p, mesh; phisolution = array_field)

    vector_fieldx, vector_fieldy, vector_fieldz = pressure_phi_gradient(velocity.p, mesh; phisolution = array_field, T = T)

    vector_to_phi!(vector_fieldx, velocity.p, mesh; phisolution = array_gradFieldx)
    vector_to_phi!(vector_fieldy, velocity.p, mesh; phisolution = array_gradFieldy)
    vector_to_phi!(vector_fieldz, velocity.p, mesh; phisolution = array_gradFieldz)

    array_div_vel = divergence_velocityToArray(
        velocity,
        mesh,
        material;
        velocityU = velocityU,
        velocityV = velocityV,
        velocityW = velocityW,
        T = T,
        interpolation = interpolation,
    )

    # Pressure
    velocity.p.eval .= array_field + pressure - array_div_vel

    # Velocity Projection
    velocity.u.eval .= velocity.u.eval - (coef * ((1.0 ./ material.ρ) .* array_gradFieldx) ./ mesh.vol)
    velocity.v.eval .= velocity.v.eval - (coef * ((1.0 ./ material.ρ) .* array_gradFieldy) ./ mesh.vol)
    velocity.w.eval .= velocity.w.eval - (coef * ((1.0 ./ material.ρ) .* array_gradFieldz) ./ mesh.vol)

    return nothing
end
