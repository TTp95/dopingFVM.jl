"""

"""
function create_Phi end

function create_Phi(
    mesh::UnionCSMesh1D;
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
)
    #allocate variables
    eval = zeros(T, mesh.l1)
    iter = zeros(T, mesh.l1)
    time1 = zeros(T, mesh.l1)
    time2 = zeros(T, mesh.l1)
    time3 = zeros(T, mesh.l1)
    sourceC = zeros(T, mesh.l1)
    sourceP = zeros(T, mesh.l1)
    onoff = ones(Bool, mesh.l1)
    gIndex = zeros(N, mesh.l1)
    bounds = zeros(Bool, mesh.l1)
    nbounds = zeros(N, mesh.l1)
    convergence = 0.0
    convergenceRelative = 0.0
    convergenceTime = 0.0
    tolerance = 0.0
    key = ""

    return CSPhi1D{T,N}(
        eval,
        iter,
        time1,
        time2,
        time3,
        sourceC,
        sourceP,
        onoff,
        gIndex,
        bounds,
        nbounds,
        convergence,
        convergenceRelative,
        convergenceTime,
        tolerance,
        key,
    )
end

function create_Phi(
    mesh::UnionCSMesh2D;
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
)
    #allocate variables
    eval = zeros(T, mesh.l1, mesh.m1)
    iter = zeros(T, mesh.l1, mesh.m1)
    time1 = zeros(T, mesh.l1, mesh.m1)
    time2 = zeros(T, mesh.l1, mesh.m1)
    time3 = zeros(T, mesh.l1, mesh.m1)
    sourceC = zeros(T, mesh.l1, mesh.m1)
    sourceP = zeros(T, mesh.l1, mesh.m1)
    onoff = ones(Bool, mesh.l1, mesh.m1)
    gIndex = zeros(N, mesh.l1, mesh.m1)
    bounds = zeros(Bool, mesh.l1, mesh.m1)
    nbounds = zeros(N, mesh.l1, mesh.m1)
    convergence = 0.0
    convergenceRelative = 0.0
    convergenceTime = 0.0
    tolerance = 0.0
    key = ""

    return CSPhi2D{T,N}(
        eval,
        iter,
        time1,
        time2,
        time3,
        sourceC,
        sourceP,
        onoff,
        gIndex,
        bounds,
        nbounds,
        convergence,
        convergenceRelative,
        convergenceTime,
        tolerance,
        key,
    )
end

function create_Phi(
    mesh::UnionCSMesh3D;
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
)
    #allocate variables
    eval = zeros(T, mesh.l1, mesh.m1, mesh.n1)
    iter = zeros(T, mesh.l1, mesh.m1, mesh.n1)
    time1 = zeros(T, mesh.l1, mesh.m1, mesh.n1)
    time2 = zeros(T, mesh.l1, mesh.m1, mesh.n1)
    time3 = zeros(T, mesh.l1, mesh.m1, mesh.n1)
    sourceC = zeros(T, mesh.l1, mesh.m1, mesh.n1)
    sourceP = zeros(T, mesh.l1, mesh.m1, mesh.n1)
    onoff = ones(Bool, mesh.l1, mesh.m1, mesh.n1)
    gIndex = zeros(N, mesh.l1, mesh.m1, mesh.n1)
    bounds = zeros(Bool, mesh.l1, mesh.m1, mesh.n1)
    nbounds = zeros(N, mesh.l1, mesh.m1, mesh.n1)
    convergence = 0.0
    convergenceRelative = 0.0
    convergenceTime = 0.0
    tolerance = 0.0
    key = ""

    return CSPhi3D{T,N}(
        eval,
        iter,
        time1,
        time2,
        time3,
        sourceC,
        sourceP,
        onoff,
        gIndex,
        bounds,
        nbounds,
        convergence,
        convergenceRelative,
        convergenceTime,
        tolerance,
        key,
    )
end

"""

"""
function create_FaceVelocity end

function create_FaceVelocity(
    mesh::UnionCSMesh1D;
    T::Type{<:AbstractFloat} = Float64,
)
    uFace = zeros(T, mesh.l1 + 1)
    uFaceIter = zeros(T, mesh.l1 + 1)
    uFaceTime1 = zeros(T, mesh.l1 + 1)
    uFaceTime2 = zeros(T, mesh.l1 + 1)
    uFaceTime3 = zeros(T, mesh.l1 + 1)

    return CSFaceVelocity1D{T}(
        uFace,
        uFaceIter,
        uFaceTime1,
        uFaceTime2,
        uFaceTime3,
    )
end

function create_FaceVelocity(
    mesh::UnionCSMesh2D;
    T::Type{<:AbstractFloat} = Float64,
)
    uFace = zeros(T, mesh.l1 + 1, mesh.m1)
    vFace = zeros(T, mesh.l1, mesh.m1 + 1)
    uFaceIter = zeros(T, mesh.l1 + 1, mesh.m1)
    vFaceIter = zeros(T, mesh.l1, mesh.m1 + 1)
    uFaceTime1 = zeros(T, mesh.l1 + 1, mesh.m1)
    vFaceTime1 = zeros(T, mesh.l1, mesh.m1 + 1)
    uFaceTime2 = zeros(T, mesh.l1 + 1, mesh.m1)
    vFaceTime2 = zeros(T, mesh.l1, mesh.m1 + 1)
    uFaceTime3 = zeros(T, mesh.l1 + 1, mesh.m1)
    vFaceTime3 = zeros(T, mesh.l1 , mesh.m1 + 1)

    return CSFaceVelocity2D{T}(
        uFace,
        vFace,
        uFaceIter,
        vFaceIter,
        uFaceTime1,
        vFaceTime1,
        uFaceTime2,
        vFaceTime2,
        uFaceTime3,
        vFaceTime3,
    )
end

function create_FaceVelocity(
    mesh::UnionCSMesh3D;
    T::Type{<:AbstractFloat} = Float64,
)
    uFace = zeros(T, mesh.l1 + 1, mesh.m1, mesh.n1)
    vFace = zeros(T, mesh.l1, mesh.m1 + 1, mesh.n1)
    wFace = zeros(T, mesh.l1, mesh.m1, mesh.n1 + 1)
    uFaceIter = zeros(T, mesh.l1 + 1, mesh.m1, mesh.n1)
    vFaceIter = zeros(T, mesh.l1, mesh.m1 + 1, mesh.n1)
    wFaceIter = zeros(T, mesh.l1, mesh.m1, mesh.n1 + 1)
    uFaceTime1 = zeros(T, mesh.l1 + 1, mesh.m1, mesh.n1)
    vFaceTime1 = zeros(T, mesh.l1, mesh.m1 + 1, mesh.n1)
    wFaceTime1 = zeros(T, mesh.l1, mesh.m1, mesh.n1 + 1)
    uFaceTime2 = zeros(T, mesh.l1 + 1, mesh.m1, mesh.n1)
    vFaceTime2 = zeros(T, mesh.l1, mesh.m1 + 1, mesh.n1)
    wFaceTime2 = zeros(T, mesh.l1, mesh.m1, mesh.n1 + 1)
    uFaceTime3 = zeros(T, mesh.l1 + 1, mesh.m1, mesh.n1)
    vFaceTime3 = zeros(T, mesh.l1, mesh.m1 + 1, mesh.n1)
    wFaceTime3 = zeros(T, mesh.l1, mesh.m1, mesh.n1 + 1)

    return CSFaceVelocity3D{T}(
        uFace,
        vFace,
        wFace,
        uFaceIter,
        vFaceIter,
        wFaceIter,
        uFaceTime1,
        vFaceTime1,
        wFaceTime1,
        uFaceTime2,
        vFaceTime2,
        wFaceTime2,
        uFaceTime3,
        vFaceTime3,
        wFaceTime3,
    )
end

"""

"""
function create_Velocity end

function create_Velocity(
    mesh::UnionCSMesh1D;
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
    template::Bool = true,
)
    #allocate variables
    u = create_Phi(mesh; T=T, N=N)
    p = create_Phi(mesh; T=T, N=N)
    fValues = create_FaceVelocity(mesh; T=T)

    return CSVelocity1D{T}(
        u,
        p,
        fValues,
    )
end

function create_Velocity(
    mesh::UnionCSMesh2D;
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
    template::Bool = true,
)
    #allocate variables
    u = create_Phi(mesh; T=T, N=N)
    v = create_Phi(mesh; T=T, N=N)
    p = create_Phi(mesh; T=T, N=N)
    fValues = create_FaceVelocity(mesh; T=T)

    return CSVelocity2D{T}(
        u,
        v,
        p,
        fValues,
    )
end

function create_Velocity(
    mesh::UnionCSMesh3D;
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
    template::Bool = true,
)
    #allocate variables
    u = create_Phi(mesh; T=T, N=N)
    v = create_Phi(mesh; T=T, N=N)
    w = create_Phi(mesh; T=T, N=N)
    p = create_Phi(mesh; T=T, N=N)
    fValues = create_FaceVelocity(mesh; T=T)

    return CSVelocity3D{T}(
        u,
        v,
        w,
        p,
        fValues,
    )
end
