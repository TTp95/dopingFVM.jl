"""
mutable struct SystemControl{T <: AbstractFloat, N <: Signed} <: TypeDopingFMV

"""
mutable struct SystemControl{T <: AbstractFloat, N <: Signed} <: TypeDopingFMV
    problem::String
    resultFolder::String
    iterate::Bool
    iterations::N
    iterationsTimeStep::N
    timeSteps::N
    print::Bool
    lastprint::T
    printTimeStep::T
    printIterations::N

end
