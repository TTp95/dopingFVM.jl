"""
This module provides:
    -

The exported names are:
$(EXPORTS)
"""
module Source

using Base.Threads

using DocStringExtensions
using SparseArrays

using dopingFVM.Structures
using dopingFVM.Tools

export discretize_source
export discretize_bodyForcesRhieChow

export implicit_relaxation!
export implicit_relaxation

include("CSSource.jl")

include("CSBodyForcesRhieChow.jl")

include("CSImplicitRelaxation.jl")

end # module
