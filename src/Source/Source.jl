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

include("CSSource.jl")

include("CSBodyForcesRhieChow.jl")

end # module
