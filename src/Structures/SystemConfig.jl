"""
mutable struct SystemControl{T <: AbstractFloat, N <: Signed} <: TypeDopingFMV

"""
mutable struct SystemConfig{T <: AbstractFloat, N <: Signed} <: TypeDopingFMV
    diffusion::N
    convection::N
    transient::N
    diffusionTime::N
    convectionTime::N
    Γntp::N
    ρntp::N
    sparrays::Bool
    mthreads::Bool
end
