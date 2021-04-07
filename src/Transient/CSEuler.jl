"""

"""
function discretize_euler_time end

function discretize_euler_time(
    phi::CSPhi1D,
    mesh::UnionCSMesh1D,
    deltat::DeltaTime,
    material::CSMaterial1D;
    materialtime1::CSMaterial1D = material,
    T::Type{<:AbstractFloat} = Float64,
    threads::Bool = false,
    sparse::Bool = true,
)
    n_equations = maximum_globalIndex(phi)

    if sparse
        A = spzeros(T, n_equations, n_equations)
    elseif !sparse
        A = zeros(T, n_equations, n_equations)
    end

    b = zeros(T, n_equations)

    for i in 1:mesh.l1
        if phi.onoff[i]
            id = phi.gIndex[i]
            num = material.ρ[i] * mesh.vol[i]
            den = deltat.dt1
            A[id,id] = num / den
            num = phi.time1[i] * materialtime1.ρ[i] * mesh.vol[i]
            den = deltat.dt1
            b[id]  = num / den
        end
    end

    return A, b
end

function discretize_euler_time(
    phi::CSPhi1D,
    mesh::UnionCSMesh1D,
    deltat::DeltaTime,
    material::UnionCSConstantMaterial;
    T::Type{<:AbstractFloat} = Float64,
    threads::Bool = false,
    sparse::Bool = true,
)
    n_equations = maximum_globalIndex(phi)

    if sparse
        A = spzeros(T, n_equations, n_equations)
    elseif !sparse
        A = zeros(T, n_equations, n_equations)
    end

    b = zeros(T, n_equations)

    for i in 1:mesh.l1
        if phi.onoff[i]
            id = phi.gIndex[i]
            num = material.ρ * mesh.vol[i]
            den = deltat.dt1
            @inbounds A[id,id]  = num / den
            num = phi.time1[i] * material.ρ * mesh.vol[i]
            den = deltat.dt1
            @inbounds b[id]  = num / den
        end
    end

    return A, b
end

function discretize_euler_time(
    phi::CSPhi2D,
    mesh::UnionCSMesh2D,
    deltat::DeltaTime,
    material::CSMaterial2D;
    materialtime1::CSMaterial2D = material,
    T::Type{<:AbstractFloat} = Float64,
    threads::Bool = false,
    sparse::Bool = true,
)
    n_equations = maximum_globalIndex(phi)

    if sparse
        A = spzeros(T, n_equations, n_equations)
    elseif !sparse
        A = zeros(T, n_equations, n_equations)
    end

    b = zeros(T, n_equations)

    for i in 1:mesh.l1
        for j in 1:mesh.m1
            if phi.onoff[i,j]
                id = phi.gIndex[i,j]
                num = material.ρ[i,j] * mesh.vol[i,j]
                den = deltat.dt1
                @inbounds A[id,id]  = num / den
                num = phi.time1[i,j] * materialtime1.ρ[i,j] * mesh.vol[i,j]
                den = deltat.dt1
                @inbounds b[id]  = num / den
            end
        end
    end

    return A, b
end

function discretize_euler_time(
    phi::CSPhi2D,
    mesh::UnionCSMesh2D,
    deltat::DeltaTime,
    material::UnionCSConstantMaterial;
    T::Type{<:AbstractFloat} = Float64,
    threads::Bool = false,
    sparse::Bool = true,
)
    n_equations = maximum_globalIndex(phi)

    if sparse
        A = spzeros(T, n_equations, n_equations)
    elseif !sparse
        A = zeros(T, n_equations, n_equations)
    end

    b = zeros(T, n_equations)

    for i in 1:mesh.l1
        for j in 1:mesh.m1
            if phi.onoff[i,j]
                id = phi.gIndex[i,j]
                num = material.ρ * mesh.vol[i,j]
                den = deltat.dt1
                @inbounds A[id,id]  = num / den
                num = phi.time1[i,j] * material.ρ * mesh.vol[i,j]
                den = deltat.dt1
                @inbounds b[id]  = num / den
            end
        end
    end

    return A, b
end

function discretize_euler_time(
    phi::CSPhi3D,
    mesh::UnionCSMesh3D,
    deltat::DeltaTime,
    material::CSMaterial3D;
    materialtime1::CSMaterial3D = material,
    T::Type{<:AbstractFloat} = Float64,
    threads::Bool = false,
    sparse::Bool = true,
    scheme::Signed = 1,
)
    n_equations = maximum_globalIndex(phi)

    if sparse
        A = spzeros(T, n_equations, n_equations)
    elseif !sparse
        A = zeros(T, n_equations, n_equations)
    end

    b = zeros(T, n_equations)

    for i in 1:mesh.l1
        for j in 1:mesh.m1
            for k in 1:mesh.n1
                if phi.onoff[i,j,k]
                    id = phi.gIndex[i,j,k]
                    num = material.ρ[i,j,k] * mesh.vol[i,j,k]
                    den = deltat.dt1
                    @inbounds A[id,id]  = num / den
                    num = phi.time1[i,j,k] * materialtime1.ρ[i,j,k] * mesh.vol[i,j,k]
                    den = deltat.dt1
                    @inbounds b[id]  = num / den
                end
            end
        end
    end

    return A, b
end

function discretize_euler_time(
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
    elseif !sparse
        A = zeros(T, n_equations, n_equations)
    end

    b = zeros(T, n_equations)

    for i in 1:mesh.l1
        for j in 1:mesh.m1
            for k in 1:mesh.n1
                if phi.onoff[i,j,k]
                    id = phi.gIndex[i,j,k]
                    num = material.ρ * mesh.vol[i,j,k]
                    den = deltat.dt1
                    @inbounds A[id,id]  = num / den
                    num = phi.time1[i,j,k] * material.ρ * mesh.vol[i,j,k]
                    den = deltat.dt1
                    @inbounds b[id]  = num / den
                end
            end
        end
    end

    return A, b
end
