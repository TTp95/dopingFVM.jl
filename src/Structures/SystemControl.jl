"""
mutable struct SystemControl{T <: AbstractFloat, N <: Signed} <: TypeDopingFMV

"""
mutable struct SystemControl{T <: AbstractFloat, N <: Signed} <: TypeDopingFMV
    problem::String
    resultFolder::String
    iterate::Bool
    iterations::N
    iterationsTimeStep::N
    minIterations::N
    maxIterations::N
    timeSteps::N
    print::Bool
    printTime::T
    printConsoleIterations::N
    lastPrintTime::T
    lastPrintConsoleIterations::T
    SolvedtimeStep::Bool
    SolvedProblem::Bool
end
