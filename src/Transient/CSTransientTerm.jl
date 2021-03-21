"""
UnionCSConstantMaterial
"""
function discretize_time end

function discretize_time(
    phi::CSPhi1D,
    mesh::UnionCSMesh1D,
    deltat::DeltaTime,
    material::CSMaterial1D;
    materialtime1::CSMaterial1D = material,
    materialtime2::CSMaterial1D = material,
    materialtime3::CSMaterial1D = material,
    T::Type{<:AbstractFloat} = Float64,
    threads::Bool = false,
    sparse::Bool = true,
    scheme::Signed = 1,
)
    n_equations = maximum_globalIndex(phi)

    if sparse
        A = spzeros(T, n_equations, n_equations)
        b = spzeros(T, n_equations)
    elseif !sparse
        A = zeros(T, n_equations, n_equations)
        b = zeros(T, n_equations)
    end

    if (scheme == 1) #euler
        if threads
            Base.Threads.@threads for i in 1:mesh.l1
                if phi.onoff[i]
                    id = phi.gIndex[i]
                    num = material.ρ[i] * mesh.vol[i]
                    den = deltat[1]
                    @inbounds A[id,id] = num / den
                    num = phi.time1[i] * materialtime1.ρ[i] * mesh.vol[i]
                    den = deltat[1]
                    @inbounds b[id]  = num / den
                end
            end
        elseif !threads
            for i in 1:mesh.l1
                if phi.onoff[i]
                    id = phi.gIndex[i]
                    num = material.ρ[i] * mesh.vol[i]
                    den = deltat[1]
                    @inbounds A[id,id] = num / den
                    num = phi.time1[i] * materialtime1.ρ[i] * mesh.vol[i]
                    den = deltat[1]
                    @inbounds b[id]  = num / den
                end
            end
        end

    elseif (scheme == 2) #CN
        a1 = ( deltat[2] / (deltat[1] * (deltat[1] + deltat[2])))
        a2 = ( (deltat[2] - deltat[1]) / (deltat[1] + deltat[2]) )
        a3 = ( deltat[1] / (deltat[2] * (deltat[1] + deltat[2])))
        if threads
            Base.Threads.@threads for i in 1:mesh.l1
                if phi.onoff[i]
                    @inbounds A[id,id] = a1 * material.ρ[i] * mesh.vol[i]
                    t1 = a2 * materialtime1.ρ[i] * phi.time1[i]
                    t2 = a3 * materialtime2.ρ[i] * phi.time2[i]
                    @inbounds b[id]  = (t1 + t2) * mesh.vol[i]
                end
            end
        elseif !threads
            for i in 1:mesh.l1
                if phi.onoff[i]
                    @inbounds A[id,id] = a1 * material.ρ[i] * mesh.vol[i]
                    t1 = a2 * materialtime1.ρ[i] * phi.time1[i]
                    t2 = a3 * materialtime2.ρ[i] * phi.time2[i]
                    @inbounds b[id]  = (t1 + t2) * mesh.vol[i]
                end
            end
        end

    elseif (scheme == 3) #BDF2
        a1 = ( 1.0 / deltat[1] ) + ( 1.0 / (deltat[1] + deltat[2]))
        a2 = (( 1.0 / deltat[1] ) + ( 1.0 / deltat[2]))
        a3 = -1.0 * ( deltat[1] / (deltat[2] * (deltat[1] + deltat[2])))
        if threads
            Base.Threads.@threads for i in 1:mesh.l1
                if phi.onoff[i]
                    @inbounds A[id,id] = a1 * material.ρ[i] * mesh.vol[i]
                    t1 = a2 * materialtime1.ρ[i] * phi.time1[i]
                    t2 = a3 * materialtime2.ρ[i] * phi.time2[i]
                    @inbounds b[id]  = (t1 + t2) * mesh.vol[i]
                end
            end
        elseif !threads
            for i in 1:mesh.l1
                if phi.onoff[i]
                    @inbounds A[id,id] = a1 * material.ρ[i] * mesh.vol[i]
                    t1 = a2 * materialtime1.ρ[i] * phi.time1[i]
                    t2 = a3 * materialtime2.ρ[i] * phi.time2[i]
                    @inbounds b[id]  = (t1 + t2) * mesh.vol[i]
                end
            end
        end

    else
        error("Time scheme: $(scheme) unimplemented")

    end

    return A, b
end

function discretize_time(
    phi::CSPhi1D,
    mesh::UnionCSMesh1D,
    deltat::DeltaTime,
    material::UnionCSConstantMaterial;
    T::Type{<:AbstractFloat} = Float64,
    threads::Bool = false,
    sparse::Bool = true,
    scheme::Signed = 1,
)
    n_equations = maximum_globalIndex(phi)

    if sparse
        A = spzeros(T, n_equations, n_equations)
        b = spzeros(T, n_equations)
    elseif !sparse
        A = zeros(T, n_equations, n_equations)
        b = zeros(T, n_equations)
    end

    if (scheme == 1) #euler
        if threads
            Base.Threads.@threads for i in 1:mesh.l1
                if phi.onoff[i]
                    id = phi.gIndex[i]
                    num = material.ρ * mesh.vol[i]
                    den = deltat[1]
                    @inbounds A[id,id]  = num / den
                    num = phi.time1[i] * material.ρ * mesh.vol[i]
                    den = deltat[1]
                    @inbounds b[id]  = num / den
                end
            end
        elseif !threads
            for i in 1:mesh.l1
                if phi.onoff[i]
                    id = phi.gIndex[i]
                    num = material.ρ * mesh.vol[i]
                    den = deltat[1]
                    @inbounds A[id,id]  = num / den
                    num = phi.time1[i] * material.ρ * mesh.vol[i]
                    den = deltat[1]
                    @inbounds b[id]  = num / den
                end
            end
        end

    elseif (scheme == 2) #CN
        a1 = ( deltat[2] / (deltat[1] * (deltat[1] + deltat[2])))
        a2 = ( (deltat[2] - deltat[1]) / (deltat[1] + deltat[2]) )
        a3 = ( deltat[1] / (deltat[2] * (deltat[1] + deltat[2])))
        if threads
            Base.Threads.@threads for i in 1:mesh.l1
                if phi.onoff[i]
                    @inbounds A[id,id]  = a1 * material.ρ * mesh.vol[i]
                    t1 = a2 * material.ρ * phi.time1[i]
                    t2 = a3 * material.ρ * phi.time2[i]
                    @inbounds b[id]  = (t1 + t2) * mesh.vol[i]
                end
            end
        elseif !threads
            for i in 1:mesh.l1
                if phi.onoff[i]
                    @inbounds A[id,id]  = a1 * material.ρ * mesh.vol[i]
                    t1 = a2 * material.ρ * phi.time1[i]
                    t2 = a3 *material.ρ * phi.time2[i]
                    @inbounds b[id]  = (t1 + t2) * mesh.vol[i]
                end
            end
        end

    elseif (scheme == 3) #BDF2
        a1 = ( 1.0 / deltat[1] ) + ( 1.0 / (deltat[1] + deltat[2]))
        a2 = (( 1.0 / deltat[1] ) + ( 1.0 / deltat[2]))
        a3 = -1.0 * ( deltat[1] / (deltat[2] * (deltat[1] + deltat[2])))
        if threads
            Base.Threads.@threads for i in 1:mesh.l1
                if phi.onoff[i]
                    @inbounds A[id,id]  = a1 * material.ρ * mesh.vol[i]
                    t1 = a2 * material.ρ * phi.time1[i]
                    t2 = a3 * material.ρ * phi.time2[i]
                    @inbounds b[id]  = (t1 + t2) * mesh.vol[i]
                end
            end
        elseif !threads
            for i in 1:mesh.l1
                if phi.onoff[i]
                    @inbounds A[id,id]  = a1 * material.ρ * mesh.vol[i]
                    t1 = a2 * material.ρ * phi.time1[i]
                    t2 = a3 * material.ρ * phi.time2[i]
                    @inbounds b[id]  = (t1 + t2) * mesh.vol[i]
                end
            end
        end

    else
        error("Time scheme: $(scheme) unimplemented")

    end

    return A, b
end

function discretize_time(
    phi::CSPhi2D,
    mesh::UnionCSMesh2D,
    deltat::DeltaTime,
    material::CSMaterial2D;
    materialtime1::CSMaterial2D = material,
    materialtime2::CSMaterial2D = material,
    materialtime3::CSMaterial2D = material,
    T::Type{<:AbstractFloat} = Float64,
    threads::Bool = false,
    sparse::Bool = true,
    scheme::Signed = 1,
)
    n_equations = maximum_globalIndex(phi)

    if sparse
        A = spzeros(T, n_equations, n_equations)
        b = spzeros(T, n_equations)
    elseif !sparse
        A = zeros(T, n_equations, n_equations)
        b = zeros(T, n_equations)
    end

    if (scheme == 1) #euler
        if threads
            Base.Threads.@threads for i in 1:mesh.l1
                for j in 1:mesh.m1
                    if phi.onoff[i,j]
                        id = phi.gIndex[i,j]
                        num = material.ρ[i,j] * mesh.vol[i,j]
                        den = deltat[1]
                        @inbounds A[id,id]  = num / den
                        num = phi.time1[i,j] * materialtime1.ρ[i,j] * mesh.vol[i,j]
                        den = deltat[1]
                        @inbounds b[id]  = num / den
                    end
                end
            end
        elseif !threads
            for i in 1:mesh.l1
                if phi.onoff[i,j]
                    id = phi.gIndex[i,j]
                    num = material.ρ[i,j] * mesh.vol[i,j]
                    den = deltat[1]
                    @inbounds A[id,id]  = num / den
                    num = phi.time1[i,j] * materialtime1.ρ[i,j] * mesh.vol[i,j]
                    den = deltat[1]
                    @inbounds b[id]  = num / den
                end
            end
        end

    elseif (scheme == 2) #CN
        a1 = ( deltat[2] / (deltat[1] * (deltat[1] + deltat[2])))
        a2 = ( (deltat[2] - deltat[1]) / (deltat[1] + deltat[2]) )
        a3 = ( deltat[1] / (deltat[2] * (deltat[1] + deltat[2])))
        if threads
            Base.Threads.@threads for i in 1:mesh.l1
                for j in 1:mesh.m1
                    if phi.onoff[i,j]
                        @inbounds A[id,id]  = a1 * material.ρ[i,j] * mesh.vol[i,j]
                        t1 = a2 * materialtime1.ρ[i,j] * phi.time1[i,j]
                        t2 = a3 * materialtime2.ρ[i,j] * phi.time2[i,j]
                        @inbounds b[id]  = (t1 + t2) * mesh.vol[i,j]
                    end
                end
            end
        elseif !threads
            for i in 1:mesh.l1
                for j in 1:mesh.m1
                    if phi.onoff[i,j]
                        @inbounds A[id,id]  = a1 * material.ρ[i,j] * mesh.vol[i,j]
                        t1 = a2 * materialtime1.ρ[i,j] * phi.time1[i,j]
                        t2 = a3 * materialtime2.ρ[i,j] * phi.time2[i,j]
                        @inbounds b[id]  = (t1 + t2) * mesh.vol[i,j]
                    end
                end
            end
        end

    elseif (scheme == 3) #BDF2
        a1 = ( 1.0 / deltat[1] ) + ( 1.0 / (deltat[1] + deltat[2]))
        a2 = (( 1.0 / deltat[1] ) + ( 1.0 / deltat[2]))
        a3 = -1.0 * ( deltat[1] / (deltat[2] * (deltat[1] + deltat[2])))
        if threads
            Base.Threads.@threads for i in 1:mesh.l1
                for j in 1:mesh.m1
                    if phi.onoff[i,j]
                        @inbounds A[id,id]  = a1 * material.ρ[i,j] * mesh.vol[i,j]
                        t1 = a2 * materialtime1.ρ[i,j] * phi.time1[i,j]
                        t2 = a3 * materialtime2.ρ[i,j] * phi.time2[i,j]
                        @inbounds b[id]  = (t1 + t2) * mesh.vol[i,j]
                    end
                end
            end
        elseif !threads
            for i in 1:mesh.l1
                if phi.onoff[i,j]
                    @inbounds A[id,id]  = a1 * material.ρ[i,j] * mesh.vol[i,j]
                    t1 = a2 * materialtime1.ρ[i,j] * phi.time1[i,j]
                    t2 = a3 * materialtime2.ρ[i,j] * phi.time2[i,j]
                    @inbounds b[id]  = (t1 + t2) * mesh.vol[i,j]
                end
            end
        end

    else
        error("Time scheme: $(scheme) unimplemented")

    end

    return A, b
end

function discretize_time(
    phi::CSPhi2D,
    mesh::UnionCSMesh2D,
    deltat::DeltaTime,
    material::UnionCSConstantMaterial;
    T::Type{<:AbstractFloat} = Float64,
    threads::Bool = false,
    sparse::Bool = true,
    scheme::Signed = 1,
)
    n_equations = maximum_globalIndex(phi)

    if sparse
        A = spzeros(T, n_equations, n_equations)
        b = spzeros(T, n_equations)
    elseif !sparse
        A = zeros(T, n_equations, n_equations)
        b = zeros(T, n_equations)
    end

    if (scheme == 1) #euler
        if threads
            Base.Threads.@threads for i in 1:mesh.l1
                for j in 1:mesh.m1
                    if phi.onoff[i,j]
                        id = phi.gIndex[i,j]
                        num = material.ρ * mesh.vol[i,j]
                        den = deltat[1]
                        @inbounds A[id,id]  = num / den
                        num = phi.time1[i,j] * material.ρ * mesh.vol[i,j]
                        den = deltat[1]
                        @inbounds b[id]  = num / den
                    end
                end
            end
        elseif !threads
            for i in 1:mesh.l1
                if phi.onoff[i,j]
                    id = phi.gIndex[i,j]
                    num = material.ρ * mesh.vol[i,j]
                    den = deltat[1]
                    @inbounds A[id,id]  = num / den
                    num = phi.time1[i,j] * material.ρ * mesh.vol[i,j]
                    den = deltat[1]
                    @inbounds b[id]  = num / den
                end
            end
        end

    elseif (scheme == 2) #CN
        a1 = ( deltat[2] / (deltat[1] * (deltat[1] + deltat[2])))
        a2 = ( (deltat[2] - deltat[1]) / (deltat[1] + deltat[2]) )
        a3 = ( deltat[1] / (deltat[2] * (deltat[1] + deltat[2])))
        if threads
            Base.Threads.@threads for i in 1:mesh.l1
                for j in 1:mesh.m1
                    if phi.onoff[i,j]
                        @inbounds A[id,id]  = a1 * material.ρ * mesh.vol[i,j]
                        t1 = a2 * material.ρ * phi.time1[i,j]
                        t2 = a3 * material.ρ * phi.time2[i,j]
                        @inbounds b[id]  = (t1 + t2) * mesh.vol[i,j]
                    end
                end
            end
        elseif !threads
            for i in 1:mesh.l1
                for j in 1:mesh.m1
                    if phi.onoff[i,j]
                        @inbounds A[id,id]  = a1 * material.ρ * mesh.vol[i,j]
                        t1 = a2 * material.ρ * phi.time1[i,j]
                        t2 = a3 * material.ρ * phi.time2[i,j]
                        @inbounds b[id]  = (t1 + t2) * mesh.vol[i,j]
                    end
                end
            end
        end

    elseif (scheme == 3) #BDF2
        a1 = ( 1.0 / deltat[1] ) + ( 1.0 / (deltat[1] + deltat[2]))
        a2 = (( 1.0 / deltat[1] ) + ( 1.0 / deltat[2]))
        a3 = -1.0 * ( deltat[1] / (deltat[2] * (deltat[1] + deltat[2])))
        if threads
            Base.Threads.@threads for i in 1:mesh.l1
                for j in 1:mesh.m1
                    if phi.onoff[i,j]
                        @inbounds A[id,id]  = a1 * material.ρ * mesh.vol[i,j]
                        t1 = a2 * material.ρ * phi.time1[i,j]
                        t2 = a3 * material.ρ * phi.time2[i,j]
                        @inbounds b[id]  = (t1 + t2) * mesh.vol[i,j]
                    end
                end
            end
        elseif !threads
            for i in 1:mesh.l1
                if phi.onoff[i,j]
                    @inbounds A[id,id]  = a1 * material.ρ * mesh.vol[i,j]
                    t1 = a2 * material.ρ * phi.time1[i,j]
                    t2 = a3 * material.ρ * phi.time2[i,j]
                    @inbounds b[id]  = (t1 + t2) * mesh.vol[i,j]
                end
            end
        end

    else
        error("Time scheme: $(scheme) unimplemented")

    end

    return A, b
end

function discretize_time(
    phi::CSPhi3D,
    mesh::UnionCSMesh3D,
    deltat::DeltaTime,
    material::CSMaterial3D;
    materialtime1::CSMaterial3D = material,
    materialtime2::CSMaterial3D = material,
    materialtime3::CSMaterial3D = material,
    T::Type{<:AbstractFloat} = Float64,
    threads::Bool = false,
    sparse::Bool = true,
    scheme::Signed = 1,
)
    n_equations = maximum_globalIndex(phi)

    if sparse
        A = spzeros(T, n_equations, n_equations)
        b = spzeros(T, n_equations)
    elseif !sparse
        A = zeros(T, n_equations, n_equations)
        b = zeros(T, n_equations)
    end

    if (scheme == 1) #euler
        if threads
            Base.Threads.@threads for i in 1:mesh.l1
                for j in 1:mesh.m1
                    for k in 1:mesh.n1
                        if phi.onoff[i,j,k]
                            id = phi.gIndex[i,j,k]
                            num = material.ρ[i,j,k] * mesh.vol[i,j,k]
                            den = deltat[1]
                            @inbounds A[id,id]  = num / den
                            num = phi.time1[i,j,k] * materialtime1.ρ[i,j,k] * mesh.vol[i,j,k]
                            den = deltat[1]
                            @inbounds b[id]  = num / den
                        end
                    end
                end
            end
        elseif !threads
            for i in 1:mesh.l1
                for j in 1:mesh.m1
                    for k in 1:mesh.n1
                        if phi.onoff[i,j,k]
                            id = phi.gIndex[i,j,k]
                            num = material.ρ[i,j,k] * mesh.vol[i,j,k]
                            den = deltat[1]
                            @inbounds A[id,id]  = num / den
                            num = phi.time1[i,j,k] * materialtime1.ρ[i,j,k] * mesh.vol[i,j,k]
                            den = deltat[1]
                            @inbounds b[id]  = num / den
                        end
                    end
                end
            end
        end

    elseif (scheme == 2) #CN
        a1 = ( deltat[2] / (deltat[1] * (deltat[1] + deltat[2])))
        a2 = ( (deltat[2] - deltat[1]) / (deltat[1] + deltat[2]) )
        a3 = ( deltat[1] / (deltat[2] * (deltat[1] + deltat[2])))
        if threads
            Base.Threads.@threads for i in 1:mesh.l1
                for j in 1:mesh.m1
                    for k in 1:mesh.n1
                        if phi.onoff[i,j,k]
                            @inbounds A[id,id]  = a1 * material.ρ[i,j,k] * mesh.vol[i,j,k]
                            t1 = a2 * materialtime1.ρ[i,j,k] * phi.time1[i,j,k]
                            t2 = a3 * materialtime2.ρ[i,j,k] * phi.time2[i,j,k]
                            @inbounds b[id]  = (t1 + t2) * mesh.vol[i,j,k]
                        end
                    end
                end
            end
        elseif !threads
            for i in 1:mesh.l1
                for j in 1:mesh.m1
                    for k in 1:mesh.n1
                        if phi.onoff[i,j,k]
                            @inbounds A[id,id]  = a1 * material.ρ[i,j,k] * mesh.vol[i,j,k]
                            t1 = a2 * materialtime1.ρ[i,j,k] * phi.time1[i,j,k]
                            t2 = a3 * materialtime2.ρ[i,j,k] * phi.time2[i,j,k]
                            @inbounds b[id]  = (t1 + t2) * mesh.vol[i,j,k]
                        end
                    end
                end
            end
        end

    elseif (scheme == 3) #BDF2
        a1 = ( 1.0 / deltat[1] ) + ( 1.0 / (deltat[1] + deltat[2]))
        a2 = (( 1.0 / deltat[1] ) + ( 1.0 / deltat[2]))
        a3 = -1.0 * ( deltat[1] / (deltat[2] * (deltat[1] + deltat[2])))
        if threads
            Base.Threads.@threads for i in 1:mesh.l1
                for j in 1:mesh.m1
                    for k in 1:mesh.n1
                        if phi.onoff[i,j,k]
                            @inbounds A[id,id]  = a1 * material.ρ[i,j,k] * mesh.vol[i,j,k]
                            t1 = a2 * materialtime1.ρ[i,j,k] * phi.time1[i,j,k]
                            t2 = a3 * materialtime2.ρ[i,j,k] * phi.time2[i,j,k]
                            @inbounds b[id]  = (t1 + t2) * mesh.vol[i,j,k]
                        end
                    end
                end
            end
        elseif !threads
            for i in 1:mesh.l1
                for j in 1:mesh.m1
                    for k in 1:mesh.n1
                        if phi.onoff[i,j,k]
                            @inbounds A[id,id]  = a1 * material.ρ[i,j,k] * mesh.vol[i,j,k]
                            t1 = a2 * materialtime1.ρ[i,j,k] * phi.time1[i,j,k]
                            t2 = a3 * materialtime2.ρ[i,j,k] * phi.time2[i,j,k]
                            @inbounds b[id]  = (t1 + t2) * mesh.vol[i,j,k]
                        end
                    end
                end
            end
        end

    else
        error("Time scheme: $(scheme) unimplemented")

    end

    return A, b
end

function discretize_time(
    phi::CSPhi3D,
    mesh::UnionCSMesh3D,
    deltat::DeltaTime,
    material::UnionCSConstantMaterial;
    T::Type{<:AbstractFloat} = Float64,
    threads::Bool = false,
    sparse::Bool = true,
    scheme::Signed = 1,
)
    n_equations = maximum_globalIndex(phi)

    if sparse
        A = spzeros(T, n_equations, n_equations)
        b = spzeros(T, n_equations)
    elseif !sparse
        A = zeros(T, n_equations, n_equations)
        b = zeros(T, n_equations)
    end

    if (scheme == 1) #euler
        if threads
            Base.Threads.@threads for i in 1:mesh.l1
                for j in 1:mesh.m1
                    for k in 1:mesh.n1
                        if phi.onoff[i,j,k]
                            id = phi.gIndex[i,j,k]
                            num = material.ρ * mesh.vol[i,j,k]
                            den = deltat[1]
                            @inbounds A[id,id]  = num / den
                            num = phi.time1[i,j,k] * material.ρ * mesh.vol[i,j,k]
                            den = deltat[1]
                            @inbounds b[id]  = num / den
                        end
                    end
                end
            end
        elseif !threads
            for i in 1:mesh.l1
                for j in 1:mesh.m1
                    for k in 1:mesh.n1
                        if phi.onoff[i,j,k]
                            id = phi.gIndex[i,j,k]
                            num = material.ρ * mesh.vol[i,j,k]
                            den = deltat[1]
                            @inbounds A[id,id]  = num / den
                            num = phi.time1[i,j,k] * material.ρ * mesh.vol[i,j,k]
                            den = deltat[1]
                            @inbounds b[id]  = num / den
                        end
                    end
                end
            end
        end

    elseif (scheme == 2) #CN
        a1 = ( deltat[2] / (deltat[1] * (deltat[1] + deltat[2])))
        a2 = ( (deltat[2] - deltat[1]) / (deltat[1] + deltat[2]) )
        a3 = ( deltat[1] / (deltat[2] * (deltat[1] + deltat[2])))
        if threads
            Base.Threads.@threads for i in 1:mesh.l1
                for j in 1:mesh.m1
                    for k in 1:mesh.n1
                        if phi.onoff[i,j,k]
                            @inbounds A[id,id]  = a1 * material.ρ * mesh.vol[i,j,k]
                            t1 = a2 * material.ρ * phi.time1[i,j,k]
                            t2 = a3 * material.ρ * phi.time2[i,j,k]
                            @inbounds b[id]  = (t1 + t2) * mesh.vol[i,j,k]
                        end
                    end
                end
            end
        elseif !threads
            for i in 1:mesh.l1
                for j in 1:mesh.m1
                    for k in 1:mesh.n1
                        if phi.onoff[i,j,k]
                            @inbounds A[id,id]  = a1 * material.ρ * mesh.vol[i,j,k]
                            t1 = a2 * material.ρ * phi.time1[i,j,k]
                            t2 = a3 * material.ρ * phi.time2[i,j,k]
                            @inbounds b[id]  = (t1 + t2) * mesh.vol[i,j,k]
                        end
                    end
                end
            end
        end

    elseif (scheme == 3) #BDF2
        a1 = ( 1.0 / deltat[1] ) + ( 1.0 / (deltat[1] + deltat[2]))
        a2 = (( 1.0 / deltat[1] ) + ( 1.0 / deltat[2]))
        a3 = -1.0 * ( deltat[1] / (deltat[2] * (deltat[1] + deltat[2])))
        if threads
            Base.Threads.@threads for i in 1:mesh.l1
                for j in 1:mesh.m1
                    for k in 1:mesh.n1
                        if phi.onoff[i,j,k]
                            @inbounds A[id,id]  = a1 * material.ρ * mesh.vol[i,j,k]
                            t1 = a2 * material.ρ * phi.time1[i,j,k]
                            t2 = a3 * material.ρ * phi.time2[i,j,k]
                            @inbounds b[id]  = (t1 + t2) * mesh.vol[i,j,k]
                        end
                    end
                end
            end
        elseif !threads
            for i in 1:mesh.l1
                for j in 1:mesh.m1
                    for k in 1:mesh.n1
                        if phi.onoff[i,j,k]
                            @inbounds A[id,id]  = a1 * material.ρ * mesh.vol[i,j,k]
                            t1 = a2 * material.ρ * phi.time1[i,j,k]
                            t2 = a3 * material.ρ * phi.time2[i,j,k]
                            @inbounds b[id]  = (t1 + t2) * mesh.vol[i,j,k]
                        end
                    end
                end
            end
        end

    else
        error("Time scheme: $(scheme) unimplemented")

    end

    return A, b
end
