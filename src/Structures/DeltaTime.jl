"""
CSPhi1D{T<:AbstractFloat, N <: Signed} <: PhiCartesianStructured

"""
mutable struct DeltaTime{T<:AbstractFloat} <: TypeDopingFMV
    Transient::Bool
    initialTime::T
    finalTime::T
    dt1::T
    dt2::T
    dt3::T
end
