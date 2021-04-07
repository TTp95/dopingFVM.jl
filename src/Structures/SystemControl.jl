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
    time::T
    timeSteps::N
    plot::Bool
    print::Bool
    inspect::Bool
    plotTime::T
    printIterations::N
    printInspect::N
    lastPlot::T
    lastPrint::N
    lastInspect::N
    solvedTimeStep::Bool
    solvedProblem::Bool
end
