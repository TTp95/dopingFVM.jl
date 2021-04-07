"""
CSPhi1D{T<:AbstractFloat, N <: Signed} <: PhiCartesianStructured

"""
mutable struct CSPhi1D{T <: AbstractFloat, N <: Signed} <: PhiCartesianStructured
    eval::Array{T,1}
    iter::Array{T,1}
    time1::Array{T,1}
    time2::Array{T,1}
    time3::Array{T,1}
    sourceC::Array{T,1}
    sourceP::Array{T,1}
    onoff::Array{Bool,1}
    gIndex::Array{<:N,1}
    bounds::Array{Bool,1}
    nbounds::Array{<:N,1}
    convergence::T
    convergenceRelative::T
    convergenceTime::T
    tolerance::T
    errorI::T
    errorR::T
    errorT::T
    residual::T
    key::String
end

"""
CSPhi2D{T<:AbstractFloat, N <: Signed} <: PhiCartesianStructured

"""
mutable struct CSPhi2D{T<:AbstractFloat, N <: Signed} <: PhiCartesianStructured
    eval::Array{T,2}
    iter::Array{T,2}
    time1::Array{T,2}
    time2::Array{T,2}
    time3::Array{T,2}
    sourceC::Array{T,2}
    sourceP::Array{T,2}
    onoff::Array{Bool,2}
    gIndex::Array{<:N,2}
    bounds::Array{Bool,2}
    nbounds::Array{<:N,2}
    convergence::T
    convergenceRelative::T
    convergenceTime::T
    tolerance::T
    errorI::T
    errorR::T
    errorT::T
    residual::T
    key::String
end

"""
CSPhi3D{T<:AbstractFloat, N <: Signed} <: PhiCartesianStructured

"""
mutable struct CSPhi3D{T<:AbstractFloat, N <: Signed} <: PhiCartesianStructured
    eval::Array{T,3}
    iter::Array{T,3}
    time1::Array{T,3}
    time2::Array{T,3}
    time3::Array{T,3}
    sourceC::Array{T,3}
    sourceP::Array{T,3}
    onoff::Array{Bool,3}
    gIndex::Array{<:N,3}
    bounds::Array{Bool,3}
    nbounds::Array{<:N,3}
    convergence::T
    convergenceRelative::T
    convergenceTime::T
    tolerance::T
    errorI::T
    errorR::T
    errorT::T
    residual::T
    key::String
end

"""
CSFaceVelocity1D{T<:AbstractFloat} <: PhiCartesianStructured

"""
struct CSFaceVelocity1D{T<:AbstractFloat} <: PhiCartesianStructured
    uFace::Array{T,1}
    uFaceIter::Array{T,1}
    uFaceTime1::Array{T,1}
    uFaceTime2::Array{T,1}
    uFaceTime3::Array{T,1}
end

"""
CSFaceVelocity2D{T<:AbstractFloat} <: PhiCartesianStructured

"""
struct CSFaceVelocity2D{T<:AbstractFloat} <: PhiCartesianStructured
    uFace::Array{T,2}
    vFace::Array{T,2}
    uFaceIter::Array{T,2}
    vFaceIter::Array{T,2}
    uFaceTime1::Array{T,2}
    vFaceTime1::Array{T,2}
    uFaceTime2::Array{T,2}
    vFaceTime2::Array{T,2}
    uFaceTime3::Array{T,2}
    vFaceTime3::Array{T,2}
end

"""
CSFaceVelocity3D{T<:AbstractFloat} <: PhiCartesianStructured

"""
struct CSFaceVelocity3D{T<:AbstractFloat} <: PhiCartesianStructured
    uFace::Array{T,3}
    vFace::Array{T,3}
    wFace::Array{T,3}
    uFaceIter::Array{T,3}
    vFaceIter::Array{T,3}
    wFaceIter::Array{T,3}
    uFaceTime1::Array{T,3}
    vFaceTime1::Array{T,3}
    wFaceTime1::Array{T,3}
    uFaceTime2::Array{T,3}
    vFaceTime2::Array{T,3}
    wFaceTime2::Array{T,3}
    uFaceTime3::Array{T,3}
    vFaceTime3::Array{T,3}
    wFaceTime3::Array{T,3}
end

"""
CSVelocity1D{T<:AbstractFloat} <: PhiCartesianStructured

"""
struct CSVelocity1D{T<:AbstractFloat} <: PhiCartesianStructured
    u::CSPhi1D
    p::CSPhi1D
    fValues::CSFaceVelocity1D
end

"""
CSVelocity2D{T<:AbstractFloat} <: PhiCartesianStructured

"""
struct CSVelocity2D{T<:AbstractFloat} <: PhiCartesianStructured
    u::CSPhi2D
    v::CSPhi2D
    p::CSPhi2D
    fValues::CSFaceVelocity2D
end

"""
CSVelocity3D{T<:AbstractFloat} <: PhiCartesianStructured

"""
struct CSVelocity3D{T<:AbstractFloat} <: PhiCartesianStructured
    u::CSPhi3D
    v::CSPhi3D
    w::CSPhi3D
    p::CSPhi3D
    fValues::CSFaceVelocity3D
end
