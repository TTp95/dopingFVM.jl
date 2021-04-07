"""

"""
function discretize_time end

function discretize_time(
    phi::UnionCSPhi,
    mesh::UnionCSMesh,
    deltat::DeltaTime,
    material::UnionCSMaterial,
    system::SystemControl;
    materialtime1::UnionCSMaterial = material,
    materialtime2::UnionCSMaterial = material,
    materialtime3::UnionCSMaterial = material,
    T::Type{<:AbstractFloat} = Float64,
    threads::Bool = false,
    sparse::Bool = true,
    scheme::Signed = 1,
    forceScheme::Bool = false,
)
    if (scheme == 1) # euler
        A, b = discretize_euler_time(
            phi,
            mesh,
            deltat,
            material;
            materialtime1 = materialtime1,
            T = T,
            threads = threads,
            sparse = sparse,
        )

    elseif (scheme == 2) # CrankNicolson
        if (system.timeSteps >= 2) || forceScheme
            A, b = discretize_crankNicolson_time(
                phi,
                mesh,
                deltat,
                material;
                materialtime1 = materialtime1,
                materialtime2 = materialtime2,
                T = T,
                threads = threads,
                sparse = sparse,
            )

        elseif (system.timeSteps == 1)
            A, b = discretize_euler_time(
                phi,
                mesh,
                deltat,
                material;
                materialtime1 = materialtime1,
                T = T,
                threads = threads,
                sparse = sparse,
            )
        end

    elseif (scheme == 3) # BDF2
        if (system.timeSteps >= 2) || forceScheme
            A, b = discretize_BDF2_time(
                phi,
                mesh,
                deltat,
                material;
                materialtime1 = materialtime1,
                materialtime2 = materialtime2,
                T = T,
                threads = threads,
                sparse = sparse,
            )

        elseif (system.timeSteps == 1)
            A, b = discretize_euler_time(
                phi,
                mesh,
                deltat,
                material;
                materialtime1 = materialtime1,
                T = T,
                threads = threads,
                sparse = sparse,
            )
        end

    else
        error("Time scheme: $(scheme) unimplemented")
    end

    return A, b
end

function discretize_time(
    phi:: UnionCSPhi,
    mesh::UnionCSMesh,
    deltat::DeltaTime,
    material::UnionCSConstantMaterial,
    system::SystemControl;
    T::Type{<:AbstractFloat} = Float64,
    threads::Bool = false,
    sparse::Bool = true,
    scheme::Signed = 1,
    forceScheme::Bool = false,
)
    if (scheme == 1) # euler
        A, b = discretize_euler_time(
            phi,
            mesh,
            deltat,
            material;
            T = T,
            threads = threads,
            sparse = sparse,
        )

    elseif (scheme == 2) # CrankNicolson
        if (system.timeSteps >= 2) || forceScheme
            A, b = discretize_crankNicolson_time(
                phi,
                mesh,
                deltat,
                material;
                T = T,
                threads = threads,
                sparse = sparse,
            )

        elseif (system.timeSteps == 1)
            A, b = discretize_euler_time(
                phi,
                mesh,
                deltat,
                material;
                T = T,
                threads = threads,
                sparse = sparse,
            )
        end

    elseif (scheme == 3) # BDF2
        if (system.timeSteps >= 2) || forceScheme
            A, b = discretize_BDF2_time(
                phi,
                mesh,
                deltat,
                material;
                T = T,
                threads = threads,
                sparse = sparse,
            )

        elseif (system.timeSteps == 1)
            A, b = discretize_euler_time(
                phi,
                mesh,
                deltat,
                material;
                T = T,
                threads = threads,
                sparse = sparse,
            )
        end

    else
        error("Time scheme: $(scheme) unimplemented")
    end

    return A, b
end
