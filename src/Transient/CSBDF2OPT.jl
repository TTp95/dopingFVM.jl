"""

"""
function discretize_BDF2OPT_time end

function discretize_BDF2OPT_time(
    phi::UnionCSPhi,
    mesh::UnionCSMesh,
    deltat::DeltaTime,
    valuetime1::AbstractArray,
    valuetime2::AbstractArray,
    valuetime3::AbstractArray,
    gindex::AbstractArray,
    onoff::AbstractArray,
    material::AbstractArray,
    materialtime1::AbstractArray = material,
    materialtime2::AbstractArray = material,
    materialtime3::AbstractArray = material;
    T::Type{<:AbstractFloat} = Float64,
    β::Float64 = 0.48,
)
    n_equations = maximum_globalIndex(phi)

    b = zeros(T, n_equations)

    a1 = β * ((1.0/deltat.dt1) * (11.0  / 6.0)) + (1 - β) * (( 1.0 / deltat.dt1 ) + ( 1.0 / (deltat.dt1 + deltat.dt2)))
    a2 = β * ((1.0/deltat.dt1) * (18.0  / 6.0)) + (1 - β) * ((( 1.0 / deltat.dt1 ) + ( 1.0 / deltat.dt2)))
    a3 = β * ((1.0/deltat.dt1) * (-9.0  / 6.0)) + (1 - β) * (-1.0 * ( deltat.dt1 / (deltat.dt2 * (deltat.dt1 + deltat.dt2))))
    a4 = β * ((1.0/deltat.dt1) * (2.0  / 6.0))

    n = 0
    D = zeros(T, n_equations)
    for i in eachindex(gindex)
        if onoff[i]
            id = gindex[i]
            n += 1
            D[n] = a1 * material[i] * mesh.vol[i]
            t1 = a2 * materialtime1[i] * valuetime1[i]
            t2 = a3 * materialtime2[i] * valuetime2[i]
            t3 = a4 * materialtime3[i] * valuetime3[i]
            b[id]  = (t1 + t2 + t3) * mesh.vol[i]
        end
    end

    A = spdiagm(0 => D[1:n])

    return A, b
end

"""

"""
function discretize_BDF2OPT_nonUniform_time end

function discretize_BDF2OPT_nonUniform_time(
    phi::UnionCSPhi,
    mesh::UnionCSMesh,
    deltat::DeltaTime,
    valuetime1::AbstractArray,
    valuetime2::AbstractArray,
    valuetime3::AbstractArray,
    gindex::AbstractArray,
    onoff::AbstractArray,
    material::AbstractArray,
    materialtime1::AbstractArray = material,
    materialtime2::AbstractArray = material,
    materialtime3::AbstractArray = material;
    T::Type{<:AbstractFloat} = Float64,
    β::Float64 = 0.48,
)
    n_equations = maximum_globalIndex(phi)

    b = zeros(T, n_equations)

    ceficents = BDF3nonUniformCoeficients(deltat)

    a1 = β * ceficents[1] + (1 - β) * (( 1.0 / deltat.dt1 ) + ( 1.0 / (deltat.dt1 + deltat.dt2)))
    a2 = β * ceficents[2] + (1 - β) * ((( 1.0 / deltat.dt1 ) + ( 1.0 / deltat.dt2)))
    a3 = β * ceficents[3] + (1 - β) * (-1.0 * ( deltat.dt1 / (deltat.dt2 * (deltat.dt1 + deltat.dt2))))
    a4 = β * ceficents[4]

    n = 0
    D = zeros(T, n_equations)
    for i in eachindex(gindex)
        if onoff[i]
            id = gindex[i]
            n += 1
            D[n] = a1 * material[i] * mesh.vol[i]
            t1 = a2 * materialtime1[i] * valuetime1[i]
            t2 = a3 * materialtime2[i] * valuetime2[i]
            t3 = a4 * materialtime3[i] * valuetime3[i]
            b[id]  = (t1 + t2 + t3) * mesh.vol[i]
        end
    end
    A = spdiagm(0 => D[1:n])

    return A, b
end
