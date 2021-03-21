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

export _globalIndex_
export assign_globalIndex!
export maximum_globalIndex

export order_iter!
export order_time!

export phi_to_vector
export vector_to_phi!
export remplace_offSolution!

export gamma_interpolation
export density_interpolation

include("CSGlobalIndex.jl")

include("CSOrderIter.jl")

include("CSOrderTime.jl")

include("CSDensityInterpolation.jl")

include("CSGammaInterpolation.jl")

include("CSPhiToVector.jl")

include("CSVectorToPhi.jl")

include("CSReemplaceOffSolution.jl")

end # module
