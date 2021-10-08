"""

"""
function discretize_BDF2_time end

function discretize_BDF2_time(
    phi::UnionCSPhi,
    mesh::UnionCSMesh,
    deltat::DeltaTime,
    valuetime1::AbstractArray,
    valuetime2::AbstractArray,
    gindex::AbstractArray,
    onoff::AbstractArray,
    material::AbstractArray,
    materialtime1::AbstractArray = material,
    materialtime2::AbstractArray = material;
    T::Type{<:AbstractFloat} = Float64,
)
    n_equations = maximum_globalIndex(phi)

    b = zeros(T, n_equations)

    a1 = ( 1.0 / deltat.dt1 ) + ( 1.0 / (deltat.dt1 + deltat.dt2))
    a2 = (( 1.0 / deltat.dt1 ) + ( 1.0 / deltat.dt2))
    a3 = -1.0 * ( deltat.dt1 / (deltat.dt2 * (deltat.dt1 + deltat.dt2)))

    n = 0
    D = zeros(T, n_equations)
    for i in eachindex(gindex)
        if onoff[i]
            id = gindex[i]
            n += 1
            D[n] = a1 * material[i] * mesh.vol[i]
            t1 = a2 * materialtime1[i] * valuetime1[i]
            t2 = a3 * materialtime2[i] * valuetime2[i]
            b[id]  = (t1 + t2) * mesh.vol[i]
        end
    end
    A = spdiagm(0 => D[1:n])

    return A, b
end
