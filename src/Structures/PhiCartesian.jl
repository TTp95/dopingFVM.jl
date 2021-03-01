"""
CSPhi1D{T<:AbstractFloat, N <: Signed} <: PhiCartesianStructured

"""
struct CSPhi1D{T<:AbstractFloat, N <: Signed} <: PhiCartesianStructured
    eval::Array{T,1}
    iter::Array{T,1}
    time1::Array{T,1}
    time2::Array{T,1}
    time3::Array{T,1}
    sourceC::Array{T,1}
    sourceP::Array{T,1}
    nn::Array{T,1}
    gIndex::Array{N,1}
    bound::Array{N,2}
end

"""
CSPhi2D{T<:AbstractFloat, N <: Signed} <: PhiCartesianStructured

"""
struct CSPhi2D{T<:AbstractFloat, N <: Signed} <: PhiCartesianStructured
    eval::Array{T,2}
    iter::Array{T,2}
    time1::Array{T,2}
    time2::Array{T,2}
    time3::Array{T,2}
    sourceC::Array{T,2}
    sourceP::Array{T,2}
    nn::Array{T,2}
    gIndex::Array{N,2}
    bound::Array{N,3}
end

"""
CSPhi3D{T<:AbstractFloat, N <: Signed} <: PhiCartesianStructured

"""
struct CSPhi3D{T<:AbstractFloat, N <: Signed} <: PhiCartesianStructured
    eval::Array{T,3}
    iter::Array{T,3}
    time1::Array{T,3}
    time2::Array{T,3}
    time3::Array{T,3}
    sourceC::Array{T,3}
    sourceP::Array{T,3}
    nn::Array{T,3}
    gIndex::Array{N,3}
    bound::Array{N,4}
end

"""
CSVelocity1D{T<:AbstractFloat, N <: Signed} <: PhiCartesianStructured

"""
struct CSVelocity1D{T<:AbstractFloat, N <: Signed} <: PhiCartesianStructured
    u::CSPhi1D
    uFace::Array{T,1}
    uFaceTime1::Array{T,1}
    uFaceTime2::Array{T,1}
    uFaceTime3::Array{T,1}
    p::CSPhi1D
end

"""
CSVelocity2D{T<:AbstractFloat, N <: Signed} <: PhiCartesianStructured

"""
struct CSVelocity2D{T<:AbstractFloat, N <: Signed} <: PhiCartesianStructured
    u::CSPhi2D
    v::CSPhi2D
    uFace::Array{T,2}
    vFace::Array{T,2}
    uFaceTime1::Array{T,2}
    vFaceTime1::Array{T,2}
    uFaceTime2::Array{T,2}
    vFaceTime2::Array{T,2}
    uFaceTime3::Array{T,2}
    vFaceTime3::Array{T,2}
    p::CSPhi2D
end

"""
CSVelocity3D{T<:AbstractFloat, N <: Signed} <: PhiCartesianStructured

"""
struct CSVelocity3D{T<:AbstractFloat, N <: Signed} <: PhiCartesianStructured
    u::CSPhi3D
    v::CSPhi3D
    w::CSPhi3D
    uFace::Array{T,3}
    vFace::Array{T,3}
    wFace::Array{T,3}
    uFaceTime1::Array{T,3}
    vFaceTime1::Array{T,3}
    wFaceTime1::Array{T,3}
    uFaceTime2::Array{T,3}
    vFaceTime2::Array{T,3}
    wFaceTime2::Array{T,3}
    uFaceTime3::Array{T,3}
    vFaceTime3::Array{T,3}
    wFaceTime3::Array{T,3}
    p::CSPhi3D
end
