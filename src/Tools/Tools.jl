"""
This module provides:
    -

The exported names are:
$(EXPORTS)
"""
module Tools

using Base.Threads

using DocStringExtensions

using dopingFVM.Structures

export assign_globalIndexCS!
export maximum_globalIndexCS
export order_iterCS!
export order_timeCS!
export remplace_offSolutionCS!
export gamma_interpolationCS
export density_interpolationCS

include("CSGlobalIndex.jl")

include("CSOrderIter.jl")

include("CSOrderTime.jl")

include("CSReemplaceOffSolution.jl")

include("CSDensityInterpolation.jl")

include("CSGammaInterpolation.jl")

end # module
