"""
mutable struct SystemControl{T <: AbstractFloat, N <: Signed} <: TypeDopingFMV

"""
mutable struct SystemConfig{T <: AbstractFloat, N <: Signed} <: TypeDopingFMV
    diffusion::N
    convection::N
    timeIntegration::N
    Γntp::N
    ρntp::N
end
