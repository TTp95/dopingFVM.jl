"""
This module provides:
    -

The exported names are:
$(EXPORTS)
"""
module VelocityInterpolation

using Base.Threads

using DocStringExtensions
using SparseArrays

using dopingFVM.Structures
using dopingFVM.Tools
using dopingFVM.Bounds

export compute_RhieChow!
export compute_RhieChow_Relaxation
export compute_RhieChow_BodyForces
export compute_RhieChow_Time

include("CSRhieChow.jl")

include("CSRhieChowRelax.jl")

include("CSRhieChowSource.jl")

include("CSRhieChowTime.jl")

include("CSRhieChowTimeConstant.jl")

end # module
