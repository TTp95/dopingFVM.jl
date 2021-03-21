"""


"""
function create_nonuniform_Mesh(
    zoneNumbers::Union{Array{<:Signed,1}, Array{<:Signed,2}},
    zoneLengths::Union{Array{<:AbstractFloat,2}, Array{<:AbstractFloat,1}},
    zoneVolumes::Union{Array{<:Signed,2}, Array{<:Signed,1}},
    zoneFactors::Union{Array{<:AbstractFloat,2}, Array{<:AbstractFloat,1}};
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
    threads::Bool = false,
    mutable::Bool = false,
)
    dimenssion = length(zoneNumbers)

    if (dimenssion==1)
        mesh = create_nonuniform_Mesh1D(
            zoneNumbers,
            zoneLengths,
            zoneVolumes,
            zoneFactors;
            T=T, N=N, threads=threads, mutable=mutable,
        )
    elseif (dimenssion==2)
        mesh = create_nonuniform_Mesh2D(
            zoneNumbers,
            zoneLengths,
            zoneVolumes,
            zoneFactors;
            T=T, N=N, threads=threads, mutable=mutable,
        )
    elseif (dimenssion==3)
        mesh = create_nonuniform_Mesh3D(
            zoneNumbers,
            zoneLengths,
            zoneVolumes,
            zoneFactors;
            T=T, N=N, threads=threads, mutable=mutable,
        )
    else
        error("input parameters generate a $(dimenssion)D mesh")
    end

    return mesh
end

"""

"""
function create_nonuniform_Mesh1D(
    zoneNumbers::Union{Array{<:Signed,1}, Array{<:Signed,2}},
    zoneLengths::Union{Array{<:AbstractFloat,2}, Array{<:AbstractFloat,1}},
    zoneVolumes::Union{Array{<:Signed,2}, Array{<:Signed,1}},
    zoneFactors::Union{Array{<:AbstractFloat,2}, Array{<:AbstractFloat,1}};
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
    threads::Bool = false,
    mutable::Bool = false,
)
    l1 = Int(sum(zoneVolumes[1,:]))       #Last volume in x-direction

    #Inicializate variables
    x = zeros(T, l1)
    dx = zeros(T, l1)
    vol = zeros(T, l1)

    #Compute each volume longitude
    @inbounds for i in 1:l1
        if (i == 1)
            global kxn = 1
            global xn = zoneVolumes[1,kxn]
            global xnn = 0
        end

        if (zoneFactors[1,kxn] != 1.0)
            num = zoneLengths[1,kxn]*(zoneFactors[1,kxn]-1.0)
            den = zoneFactors[1,kxn]^(zoneVolumes[1,kxn])-1.0
            A = num/den
            dx[i] = A*zoneFactors[1,kxn]^(i-xnn-1)
        else
            dx[i] = zoneLengths[1,kxn]/zoneVolumes[1,kxn]
        end

        if ((i == (xn)) && (kxn < zoneNumbers[1]))
            global xnn = xnn + zoneVolumes[1,kxn]
            global kxn = kxn + 1
            global xn = xn + zoneVolumes[1,kxn]
        end
    end

    #Compute nodes position in x-direction
    for i in 1:l1
        if (i == 1)
            @inbounds global x[i] = 0.5 * dx[i]
        else
            @inbounds global x[i] = x[i-1] + 0.5 * (dx[i-1] + dx[i])
        end
    end

    #Compute the volume of each control volume
    if !threads
        for i in 1:l1
            @inbounds global vol[i] = dx[i]
        end
    elseif threads
        Base.Threads.@threads for i in 1:l1
            @inbounds global vol[i] = dx[i]
        end
    end

    if mutable
        return CSMesh1D{T,N}(x, dx, vol, l1)
    elseif !mutable
        return CSMesh1DImmutable{T,N}(x, dx, vol, l1)
    end
end

"""

"""
function create_nonuniform_Mesh2D(
    zoneNumbers::Union{Array{<:Signed,1}, Array{<:Signed,2}},
    zoneLengths::Union{Array{<:AbstractFloat,2}, Array{<:AbstractFloat,1}},
    zoneVolumes::Union{Array{<:Signed,2}, Array{<:Signed,1}},
    zoneFactors::Union{Array{<:AbstractFloat,2}, Array{<:AbstractFloat,1}};
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
    threads::Bool = false,
    mutable::Bool = false,
)
    l1 = Int(sum(zoneVolumes[1,:]))       #Last volume in x-direction
    m1 = Int(sum(zoneVolumes[2,:]))       #Last volume in y-direction

    #Inicializate variables
    x = zeros(T, l1)
    y = zeros(T, m1)
    dx = zeros(T, l1)
    dy = zeros(T, m1)
    vol = zeros(T, l1, m1)

    #Compute each volume longitude
    @inbounds for i in 1:l1
        if (i == 1)
            global kxn = 1
            global xn = zoneVolumes[1,kxn]
            global xnn = 0
        end

        if (zoneFactors[1,kxn] != 1.0)
            num = zoneLengths[1,kxn]*(zoneFactors[1,kxn]-1.0)
            den = zoneFactors[1,kxn]^(zoneVolumes[1,kxn])-1.0
            A = num/den
            dx[i] = A*zoneFactors[1,kxn]^(i-xnn-1)
        else
            dx[i] = zoneLengths[1,kxn]/zoneVolumes[1,kxn]
        end

        if ((i == (xn)) && (kxn < zoneNumbers[1]))
            global xnn = xnn + zoneVolumes[1,kxn]
            global kxn = kxn + 1
            global xn = xn + zoneVolumes[1,kxn]
        end
    end

    @inbounds for j in 1:m1
        if (j == 1)
            global kyn = 1
            global yn = zoneVolumes[2,kyn]
            global ynn = 0
        end

        if (zoneFactors[2,kyn] != 1.0)
            num = zoneLengths[2,kyn]*(zoneFactors[2,kyn]-1.0)
            den = zoneFactors[2,kyn]^(zoneVolumes[2,kyn])-1.0
            A = num/den
            dy[j] = A*zoneFactors[2,kyn]^(j-ynn-1)
        else
            dy[j] = zoneLengths[2,kyn]/zoneVolumes[2,kyn]
        end

        if ((j == (yn)) && (kyn < zoneNumbers[2]))
            global ynn = ynn + zoneVolumes[2,kyn]
            global kyn = kyn + 1
            global yn = yn + zoneVolumes[2,kyn]
        end
    end

    #Compute nodes position in x-direction
    for i in 1:l1
        if (i == 1)
            @inbounds global x[i] = 0.5 * dx[i]
        else
            @inbounds global x[i] = x[i-1] + 0.5 * (dx[i-1] + dx[i])
        end
    end

    #Compute nodes position in y-direction
    for j in 1:m1
        if (j == 1)
            @inbounds global y[j] = 0.5 * dy[j]
        else
            @inbounds global y[j] = y[j-1] + 0.5 * (dy[j-1] + dy[j])
        end
    end

    #Compute the volume of each control volume
    if !threads
        for i in 1:l1
            for j in 1:m1
                @inbounds global vol[i,j] = dx[i] * dy[j]
            end
        end
    elseif threads
        Base.Threads.@threads for i in 1:l1
            for j in 1:m1
                @inbounds global vol[i,j] = dx[i] * dy[j]
            end
        end
    end

    if mutable
        return CSMesh2D{T,N}(x, y, dx, dy, vol, l1, m1)
    elseif !mutable
        return CSMesh2DImmutable{T,N}(x, y, dx, dy, vol, l1, m1)
    end
end

"""

"""
function create_nonuniform_Mesh3D(
    zoneNumbers::Union{Array{<:Signed,1}, Array{<:Signed,2}},
    zoneLengths::Union{Array{<:AbstractFloat,2}, Array{<:AbstractFloat,1}},
    zoneVolumes::Union{Array{<:Signed,2}, Array{<:Signed,1}},
    zoneFactors::Union{Array{<:AbstractFloat,2}, Array{<:AbstractFloat,1}};
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
    threads::Bool = false,
    mutable::Bool = false,
)
    l1 = Int(sum(zoneVolumes[1,:]))       #Last volume in x-direction
    m1 = Int(sum(zoneVolumes[2,:]))       #Last volume in y-direction
    n1 = Int(sum(zoneVolumes[3,:]))       #Last volume in z-direction

    #Inicializate variables
    x = zeros(T, l1)
    y = zeros(T, m1)
    z = zeros(T, n1)
    dx = zeros(T, l1)
    dy = zeros(T, m1)
    dz = zeros(T, n1)
    vol = zeros(T, l1, m1, n1)

    #Compute each volume longitude
    @inbounds for i in 1:l1
        if (i == 1)
            global kxn = 1
            global xn = zoneVolumes[1,kxn]
            global xnn = 0
        end

        if (zoneFactors[1,kxn] != 1.0)
            num = zoneLengths[1,kxn]*(zoneFactors[1,kxn]-1.0)
            den = zoneFactors[1,kxn]^(zoneVolumes[1,kxn])-1.0
            A = num/den
            dx[i] = A*zoneFactors[1,kxn]^(i-xnn-1)
        else
            dx[i] = zoneLengths[1,kxn]/zoneVolumes[1,kxn]
        end

        if ((i == (xn)) && (kxn < zoneNumbers[1]))
            global xnn = xnn + zoneVolumes[1,kxn]
            global kxn = kxn + 1
            global xn = xn + zoneVolumes[1,kxn]
        end
    end

    @inbounds for j in 1:m1
        if (j == 1)
            global kyn = 1
            global yn = zoneVolumes[2,kyn]
            global ynn = 0
        end

        if (zoneFactors[2,kyn] != 1.0)
            num = zoneLengths[2,kyn]*(zoneFactors[2,kyn]-1.0)
            den = zoneFactors[2,kyn]^(zoneVolumes[2,kyn])-1.0
            A = num/den
            dy[j] = A*zoneFactors[2,kyn]^(j-ynn-1)
        else
            dy[j] = zoneLengths[2,kyn]/zoneVolumes[2,kyn]
        end

        if ((j == (yn)) && (kyn < zoneNumbers[2]))
            global ynn = ynn + zoneVolumes[2,kyn]
            global kyn = kyn + 1
            global yn = yn + zoneVolumes[2,kyn]
        end
    end

    @inbounds for k in 1:n1
        if (k == 1)
            global kzn = 1
            global zn = zoneVolumes[3,kzn]
            global znn = 0
        end

        if (zoneFactors[3,kzn] != 1.0)
            num = zoneLengths[3,kzn]*(zoneFactors[3,kzn]-1.0)
            den = zoneFactors[3,kzn]^(zoneVolumes[3,kzn])-1.0
            A = num/den
            dz[k] = A*zoneFactors[3,kzn]^(k-znn-1)
        else
            dz[k] = zoneLengths[3,kzn]/zoneVolumes[3,kzn]
        end

        if ((k == (zn)) && (kzn < zoneNumbers[3]))
            global znn = znn + zoneVolumes[3,kzn]
            global kzn = kzn + 1
            global zn = zn + zoneVolumes[3,kyn]
        end
    end

    #Compute nodes position in x-direction
    for i in 1:l1
        if (i == 1)
            @inbounds global x[i] = 0.5 * dx[i]
        else
            @inbounds global x[i] = x[i-1] + 0.5 * (dx[i-1] + dx[i])
        end
    end

    #Compute nodes position in y-direction
    for j in 1:m1
        if (j == 1)
            @inbounds global y[j] = 0.5 * dy[j]
        else
            @inbounds global y[j] = y[j-1] + 0.5 * (dy[j-1] + dy[j])
        end
    end

    #Compute nodes position in z-direction
    for k in 1:n1
        if (k == 1)
            @inbounds global z[k] = 0.5 * dz[k]
        else
            @inbounds global z[k] = z[k-1] + 0.5 * (dz[k-1] + dz[k])
        end
    end

    #Compute the volume of each control volume
    if !threads
        for i in 1:l1
            for j in 1:m1
                for k in 1:n1
                    @inbounds global vol[i,j,k] = dx[i] * dy[j] * dz[k]
                end
            end
        end
    elseif threads
        Base.Threads.@threads for i in 1:l1
            for j in 1:m1
                for k in 1:n1
                    @inbounds global vol[i,j,k] = dx[i] * dy[j] * dz[k]
                end
            end
        end
    end

    if mutable
        return CSMesh3D{T,N}(x, y, z, dx, dy, dz, vol, l1, m1, n1)
    elseif !mutable
        return CSMesh3DImmutable{T,N}(x, y, z, dx, dy, dz, vol, l1, m1, n1)
    end
end
