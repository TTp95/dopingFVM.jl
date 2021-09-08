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
    mthreads::Bool = false,
    sparrays::Bool = true,
    scheme::Signed = 1,
    forceScheme::Bool = false,
    nonUniform::Bool = false,
)
    if (scheme == 1) # euler
        A, b = discretize_euler_time(
            phi,
            mesh,
            deltat,
            phi.time1,
            phi.gIndex,
            phi.onoff,
            material.ρ,
            materialtime1.ρ;
            T = T,
            mthreads = mthreads,
            sparrays = sparrays,
        )

    elseif (scheme == 2) # CrankNicolson
        if (system.timeSteps >= 2) || forceScheme
            A, b = discretize_crankNicolson_time(
                phi,
                mesh,
                deltat,
                phi.time1,
                phi.time2,
                phi.gIndex,
                phi.onoff,
                material.ρ,
                materialtime1.ρ,
                materialtime2.ρ;
                T = T,
                mthreads = mthreads,
                sparrays = sparrays,
            )

        elseif (system.timeSteps == 1)
            A, b = discretize_euler_time(
                phi,
                mesh,
                deltat,
                phi.time1,
                phi.gIndex,
                phi.onoff,
                material.ρ,
                materialtime1.ρ;
                T = T,
                mthreads = mthreads,
                sparrays = sparrays,
            )
        end

    elseif (scheme == 3) # BDF2
        if (system.timeSteps >= 2) || forceScheme
            A, b = discretize_BDF2_time(
                phi,
                mesh,
                deltat,
                phi.time1,
                phi.time2,
                phi.gIndex,
                phi.onoff,
                material.ρ,
                materialtime1.ρ,
                materialtime2.ρ;
                T = T,
                mthreads = mthreads,
                sparrays = sparrays,
            )

        elseif (system.timeSteps == 1)
            A, b = discretize_euler_time(
                phi,
                mesh,
                deltat,
                phi.time1,
                phi.gIndex,
                phi.onoff,
                material.ρ,
                materialtime1.ρ;
                T = T,
                mthreads = mthreads,
                sparrays = sparrays,
            )
        end

    elseif ((scheme == 4) && (!nonUniform))# BDF3 uniform time step
        if (system.timeSteps >= 3) || forceScheme
            A, b = discretize_BDF3_time(
                phi,
                mesh,
                deltat,
                phi.time1,
                phi.time2,
                phi.time3,
                phi.gIndex,
                phi.onoff,
                material.ρ,
                materialtime1.ρ,
                materialtime2.ρ,
                materialtime3.ρ;
                T = T,
                mthreads = mthreads,
                sparrays = sparrays,
            )

        elseif (system.timeSteps == 2)
            A, b = discretize_BDF2_time(
                phi,
                mesh,
                deltat,
                phi.time1,
                phi.time2,
                phi.gIndex,
                phi.onoff,
                material.ρ,
                materialtime1.ρ,
                materialtime2.ρ;
                T = T,
                mthreads = mthreads,
                sparrays = sparrays,
            )

        elseif (system.timeSteps == 1)
            A, b = discretize_euler_time(
                phi,
                mesh,
                deltat,
                phi.time1,
                phi.gIndex,
                phi.onoff,
                material.ρ,
                materialtime1.ρ;
                T = T,
                mthreads = mthreads,
                sparrays = sparrays,
            )
        end

    elseif ((scheme == 4) && (nonUniform))# BDF3 + non uniform time step
        if (system.timeSteps >= 3) || forceScheme
            A, b = discretize_BDF3_nonUniform_time(
                phi,
                mesh,
                deltat,
                phi.time1,
                phi.time2,
                phi.time3,
                phi.gIndex,
                phi.onoff,
                material.ρ,
                materialtime1.ρ,
                materialtime2.ρ,
                materialtime3.ρ;
                T = T,
                mthreads = mthreads,
                sparrays = sparrays,
            )

        elseif (system.timeSteps == 2)
            A, b = discretize_BDF2_time(
                phi,
                mesh,
                deltat,
                phi.time1,
                phi.time2,
                phi.gIndex,
                phi.onoff,
                material.ρ,
                materialtime1.ρ,
                materialtime2.ρ;
                T = T,
                mthreads = mthreads,
                sparrays = sparrays,
            )

        elseif (system.timeSteps == 1)
            A, b = discretize_euler_time(
                phi,
                mesh,
                deltat,
                phi.time1,
                phi.gIndex,
                phi.onoff,
                material.ρ,
                materialtime1.ρ;
                T = T,
                mthreads = mthreads,
                sparrays = sparrays,
            )
        end

    else
        error("Time scheme: $(scheme) unimplemented")
    end

    return A, b
end

function discretize_time(
    phi::UnionCSPhi,
    mesh::UnionCSMesh,
    deltat::DeltaTime,
    material::UnionCSConstantMaterial,
    system::SystemControl;
    T::Type{<:AbstractFloat} = Float64,
    mthreads::Bool = false,
    sparrays::Bool = true,
    scheme::Signed = 1,
    forceScheme::Bool = false,
    nonUniform::Bool = false,
)
    if (scheme == 1) # euler
        A, b = discretize_euler_time(
            phi,
            mesh,
            deltat,
            phi.time1,
            phi.gIndex,
            phi.onoff,
            material.ρ;
            T = T,
            mthreads = mthreads,
            sparrays = sparrays,
        )

    elseif (scheme == 2) # CrankNicolson
        if (system.timeSteps >= 2) || forceScheme
            A, b = discretize_crankNicolson_time(
                phi,
                mesh,
                deltat,
                phi.time1,
                phi.time2,
                phi.gIndex,
                phi.onoff,
                material.ρ;
                T = T,
                mthreads = mthreads,
                sparrays = sparrays,
            )

        elseif (system.timeSteps == 1)
            A, b = discretize_euler_time(
                phi,
                mesh,
                deltat,
                phi.time1,
                phi.gIndex,
                phi.onoff,
                material.ρ;
                T = T,
                mthreads = mthreads,
                sparrays = sparrays,
            )
        end

    elseif (scheme == 3) # BDF2
        if (system.timeSteps >= 2) || forceScheme
            A, b = discretize_BDF2_time(
                phi,
                mesh,
                deltat,
                phi.time1,
                phi.time2,
                phi.gIndex,
                phi.onoff,
                material.ρ;
                T = T,
                mthreads = mthreads,
                sparrays = sparrays,
            )

        elseif (system.timeSteps == 1)
            A, b = discretize_euler_time(
                phi,
                mesh,
                deltat,
                phi.time1,
                phi.gIndex,
                phi.onoff,
                material.ρ;
                T = T,
                mthreads = mthreads,
                sparrays = sparrays,
            )
        end

    elseif ((scheme == 4) && (!nonUniform))# BDF3 uniform time step
        if (system.timeSteps >= 3) || forceScheme
            A, b = discretize_BDF3_time(
                phi,
                mesh,
                deltat,
                phi.time1,
                phi.time2,
                phi.time3,
                phi.gIndex,
                phi.onoff,
                material.ρ;
                T = T,
                mthreads = mthreads,
                sparrays = sparrays,
            )

        elseif (system.timeSteps == 2)
            A, b = discretize_BDF2_time(
                phi,
                mesh,
                deltat,
                phi.time1,
                phi.time2,
                phi.gIndex,
                phi.onoff,
                material.ρ;
                T = T,
                mthreads = mthreads,
                sparrays = sparrays,
            )

        elseif (system.timeSteps == 1)
            A, b = discretize_euler_time(
                phi,
                mesh,
                deltat,
                phi.time1,
                phi.gIndex,
                phi.onoff,
                material.ρ;
                T = T,
                mthreads = mthreads,
                sparrays = sparrays,
            )
        end

    elseif ((scheme == 4) && (nonUniform))# BDF3 + non uniform time step
        if (system.timeSteps >= 3) || forceScheme
            A, b = discretize_BDF3_nonUniform_time(
                phi,
                mesh,
                deltat,
                phi.time1,
                phi.time2,
                phi.time3,
                phi.gIndex,
                phi.onoff,
                material.ρ;
                T = T,
                mthreads = mthreads,
                sparrays = sparrays,
            )

        elseif (system.timeSteps == 2)
            A, b = discretize_BDF2_time(
                phi,
                mesh,
                deltat,
                phi.time1,
                phi.time2,
                phi.gIndex,
                phi.onoff,
                material.ρ;
                T = T,
                mthreads = mthreads,
                sparrays = sparrays,
            )

        elseif (system.timeSteps == 1)
            A, b = discretize_euler_time(
                phi,
                mesh,
                deltat,
                phi.time1,
                phi.gIndex,
                phi.onoff,
                material.ρ;
                T = T,
                mthreads = mthreads,
                sparrays = sparrays,
            )
        end
    else
        error("Time scheme: $(scheme) unimplemented")
    end

    return A, b
end
