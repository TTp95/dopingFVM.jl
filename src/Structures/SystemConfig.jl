"""
mutable struct SystemControl{T <: AbstractFloat, N <: Signed} <: TypeDopingFMV

"""
mutable struct SystemConfig{T <: AbstractFloat, N <: Signed} <: TypeDopingFMV
    diffusion::String
    convection::String
    timeIntegration::String
end
