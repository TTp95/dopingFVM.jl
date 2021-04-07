"""

"""
function velocityProjection_PPC_Incremental! end

function velocityProjection_PPC_Incremental!(
    field::Array{<:AbstractFloat,1},
    velocity::CSVelocity1D,
    mesh::UnionCSMesh1D,
    deltat::DeltaTime,
    material::UnionCSMaterialAll;
    T::Type{<:AbstractFloat} = Float64,
    sparrays::Bool = true,
    mthreads::Bool = false,
    transientScheme::Signed = 1,
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

    vector_to_phi!(field, velocity.p, mesh; phisolution = array_field, T = T, threads = threads)

    vector_fieldx = pressure_phi_gradient(velocity.p, mesh; phisolution = array_field, T = T, threads = mthreads)

    vector_to_phi!(vector_fieldx, velocity.p, mesh; phisolution = array_gradFieldx, T = T, threads = mthreads)

    # Pressure
    velocity.p.eval .= array_field + velocity.p.time1

    # Velocity Projection
    velocity.u.eval .= velocity.u.eval - (coef * (material.ρ .* array_gradFieldx))

    return nothing
end

function velocityProjection_PPC_Incremental!(
    field::Array{<:AbstractFloat,2},
    velocity::CSVelocity2D,
    mesh::UnionCSMesh2D,
    deltat::DeltaTime,
    material::UnionCSMaterialAll;
    T::Type{<:AbstractFloat} = Float64,
    sparrays::Bool = true,
    mthreads::Bool = false,
    transientScheme::Signed = 1,
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

    vector_to_phi!(field, velocity.p, mesh; phisolution = array_field, T = T, threads = threads)

    vector_fieldx, vector_fieldy = pressure_phi_gradient(velocity.p, mesh; phisolution = array_field, T = T, threads = mthreads)

    vector_to_phi!(vector_fieldx, velocity.p, mesh; phisolution = array_gradFieldx, T = T, threads = mthreads)
    vector_to_phi!(vector_fieldy, velocity.p, mesh; phisolution = array_gradFieldy, T = T, threads = mthreads)

    # Pressure
    velocity.p.eval .= array_field + velocity.p.time1

    # Velocity Projection
    velocity.u.eval .= velocity.u.eval - (coef * (material.ρ .* array_gradFieldx))
    velocity.v.eval .= velocity.v.eval - (coef * (material.ρ .* array_gradFieldy))

    return nothing
end

function velocityProjection_PPC_Incremental!(
    field::Array{<:AbstractFloat,3},
    velocity::CSVelocity3D,
    mesh::UnionCSMesh3D,
    deltat::DeltaTime,
    material::UnionCSMaterialAll;
    T::Type{<:AbstractFloat} = Float64,
    sparrays::Bool = true,
    mthreads::Bool = false,
    transientScheme::Signed = 1,
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

    vector_to_phi!(field, velocity.p, mesh; phisolution = array_field, T = T, threads = threads)

    vector_fieldx, vector_fieldy, vector_fieldz = pressure_phi_gradient(velocity.p, mesh; phisolution = array_field, T = T, threads = mthreads)

    vector_to_phi!(vector_fieldx, velocity.p, mesh; phisolution = array_gradFieldx, T = T, threads = mthreads)
    vector_to_phi!(vector_fieldy, velocity.p, mesh; phisolution = array_gradFieldy, T = T, threads = mthreads)
    vector_to_phi!(vector_fieldz, velocity.p, mesh; phisolution = array_gradFieldz, T = T, threads = mthreads)

    # Pressure
    velocity.p.eval .= array_field + velocity.p.time1

    # Velocity Projection
    velocity.u.eval .= velocity.u.eval - (coef * (material.ρ .* array_gradFieldx))
    velocity.v.eval .= velocity.v.eval - (coef * (material.ρ .* array_gradFieldy))
    velocity.w.eval .= velocity.w.eval - (coef * (material.ρ .* array_gradFieldz))

    return nothing
end
