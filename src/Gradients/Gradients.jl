"""
This module provides:
    -

The exported names are:
$(EXPORTS)
"""
module Gradients

using Base.Threads

using DocStringExtensions

using dopingFVM.Structures
using dopingFVM.Tools
using dopingFVM.Bounds

export pressure_phi_gradient
export array_gradient
export array_coordGradient
export array_gradient_bound

include("CSPressurePhiGradient.jl")

include("CSArrayGradient.jl")

include("CSArrayCoordGradient.jl")

end # module
