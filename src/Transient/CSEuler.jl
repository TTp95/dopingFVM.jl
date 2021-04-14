"""

"""
function discretize_euler_time end

function discretize_euler_time(
    phi::UnionCSPhi,
    mesh::UnionCSMesh,
    deltat::DeltaTime,
    valuetime1::AbstractArray,
    gindex::AbstractArray,
    onoff::AbstractArray,
    material::AbstractArray,
    materialtime1::AbstractArray = material;
    T::Type{<:AbstractFloat} = Float64,
    mthreads::Bool = false,
    sparrays::Bool = true,
)
    n_equations = maximum_globalIndex(phi)

    b = zeros(T, n_equations)

    if mthreads

    elseif sparrays
        n = 0
        D = zeros(T, n_equations)
        for i in eachindex(gindex)
            if onoff[i]
                id = gindex[i]
                n += 1
                num = material[i] * mesh.vol[i]
                den = deltat.dt1
                D[n] = num / den
                num = valuetime1[i] * materialtime1[i] * mesh.vol[i]
                den = deltat.dt1
                b[id]  = num / den
            end
        end
        A = spdiagm(0 => D[1:n])
    elseif !sparrays
        A = zeros(T, n_equations, n_equations)
        for i in eachindex(gindex)
            if onoff[i]
                id = gindex[i]
                num = material[i] * mesh.vol[i]
                den = deltat.dt1
                A[id, id] = num / den
                num = valuetime1[i] * materialtime1[i] * mesh.vol[i]
                den = deltat.dt1
                b[id]  = num / den
            end
        end
    end

    return A, b
end

function discretize_euler_time(
    phi::UnionCSPhi,
    mesh::UnionCSMesh,
    deltat::DeltaTime,
    valuetime1::AbstractArray,
    gindex::AbstractArray,
    onoff::AbstractArray,
    material::AbstractFloat;
    T::Type{<:AbstractFloat} = Float64,
    mthreads::Bool = false,
    sparrays::Bool = true,
)
    n_equations = maximum_globalIndex(phi)

    b = zeros(T, n_equations)

    if mthreads

    elseif sparrays
        n = 0
        D = zeros(T, n_equations)
        for i in eachindex(gindex)
            if onoff[i]
                id = gindex[i]
                n += 1
                num = material * mesh.vol[i]
                den = deltat.dt1
                D[n] = num / den
                num = valuetime1[i] * material * mesh.vol[i]
                den = deltat.dt1
                b[id]  = num / den
            end
        end
        A = spdiagm(0 => D[1:n])
    elseif !sparrays
        A = zeros(T, n_equations, n_equations)
        for i in eachindex(gindex)
            if onoff[i]
                id = gindex[i]
                num = material * mesh.vol[i]
                den = deltat.dt1
                A[id, id] = num / den
                num = valuetime1[i] * material * mesh.vol[i]
                den = deltat.dt1
                b[id]  = num / den
            end
        end
    end

    return A, b
end
