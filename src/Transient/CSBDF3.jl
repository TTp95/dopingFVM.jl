"""

"""
function discretize_BDF3_time end

function discretize_BDF3_time(
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
)
    n_equations = maximum_globalIndex(phi)

    b = zeros(T, n_equations)

    a1 = (1.0/deltat.dt1) * (11.0  / 6.0)
    a2 = (1.0/deltat.dt1) * (18.0  / 6.0)
    a3 = (1.0/deltat.dt1) * (-9.0  / 6.0)
    a4 = (1.0/deltat.dt1) * (2.0  / 6.0)

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

"""

"""
function BDF3nonUniformCoeficients(
    deltat::DeltaTime,
)
    a = deltat.dt1
    b = a + deltat.dt2
    c = b + deltat.dt3

    AA = [1.0 1.0  1.0 1.0;
          0.0 -a -b -c;
          0.0 0.5(-a)^2 0.5(-b)^2 0.5(-c)^2;
          0.0 (1.0/6.0)*(-a)^3 (1.0/6.0)*(-b)^3 (1.0/6.0)*(-c)^3;]

    bb = [0.0; 1.0; 0.0; 0.0]

    result = (AA\bb) .* [1.0; -1.0; -1.0; -1.0]

    return result
end


"""

"""
function discretize_BDF3_nonUniform_time end

function discretize_BDF3_nonUniform_time(
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
)
    n_equations = maximum_globalIndex(phi)

    b = zeros(T, n_equations)

    ceficents = BDF3nonUniformCoeficients(deltat)

    a1 = ceficents[1]
    a2 = ceficents[2]
    a3 = ceficents[3]
    a4 = ceficents[4]

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
