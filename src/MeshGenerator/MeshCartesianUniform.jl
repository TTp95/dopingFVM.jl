"""

"""
function create_uniform_Mesh(
    lengths::Array{<:AbstractFloat,1},
    volumes::Array{<:Signed,1};
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
    threads::Bool = false,
    mutable::Bool = false,
)
    dimenssion = length(lengths)

    if (dimenssion==1)
        mesh = create_uniform_Mesh1D(
            lengths, volumes;
            T=T, N=N, threads=threads, mutable=mutable,
        )
    elseif (dimenssion==2)
        mesh = create_uniform_Mesh2D(
            lengths, volumes;
            T=T, N=N, threads=threads, mutable=mutable,
        )
    elseif (dimenssion==3)
        mesh = create_uniform_Mesh3D(
            lengths, volumes;
            T=T, N=N, threads=threads, mutable=mutable,
        )
    else
        error("input parameters generate a $(dimenssion)D mesh")
    end

    return mesh
end

"""

"""
function create_uniform_Mesh1D(
    lengths::Array{<:AbstractFloat,1},
    volumes::Array{<:Signed,1};
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
    threads::Bool = false,
    mutable::Bool = false,
)
    #Inicializate variables
    x = zeros(T, volumes[1])
    dx = zeros(T, volumes[1])
    vol = zeros(T, volumes[1])

    #Last volume in x-direction
    l1 = volumes[1]

    #Compute each volume longitude
    for i in 1:l1
        @inbounds global dx[i] = lengths[1]/volumes[1]
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
function create_uniform_Mesh2D(
    lengths::Array{<:AbstractFloat,1},
    volumes::Array{<:Signed,1};
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
    threads::Bool = false,
    mutable::Bool = false,
)
    #Inicializate variables
    x = zeros(T, volumes[1])
    y = zeros(T, volumes[2])
    dx = zeros(T, volumes[1])
    dy = zeros(T, volumes[2])
    vol = zeros(T, volumes[1], volumes[2])

    #Last volume in x-direction
    l1 = volumes[1]

    #Last volume in y-direction
    m1 = volumes[2]

    #Compute each volume longitude
    for i in 1:l1
        @inbounds global dx[i] = lengths[1]/volumes[1]
    end

    for j in 1:m1
        @inbounds global dy[j] = lengths[2]/volumes[2]
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
function create_uniform_Mesh3D(
    lengths::Array{<:AbstractFloat,1},
    volumes::Array{<:Signed,1};
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
    threads::Bool = false,
    mutable::Bool = false,
)
    #Inicializate variables
    x = zeros(T, volumes[1])
    y = zeros(T, volumes[2])
    z = zeros(T, volumes[3])
    dx = zeros(T, volumes[1])
    dy = zeros(T, volumes[2])
    dz = zeros(T, volumes[3])
    vol = zeros(T, volumes[1], volumes[2], volumes[3])

    #Last volume in x-direction
    l1 = volumes[1]

    #Last volume in y-direction
    m1 = volumes[2]

    #Last volume in z-direction
    n1 = volumes[3]

    #Compute each volume longitude
    for i in 1:l1
        @inbounds global dx[i] = lengths[1]/volumes[1]
    end

    for j in 1:m1
        @inbounds global dy[j] = lengths[2]/volumes[2]
    end

    for k in 1:n1
        @inbounds global dz[k] = lengths[3]/volumes[3]
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
