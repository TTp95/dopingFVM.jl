"""
This module provides:
    -

The exported names are:
$(EXPORTS)
"""
module Convergence

using Base.Threads

using DocStringExtensions

using dopingFVM.Structures
using dopingFVM.Tools

export convergence_iter
export convergence_relative_iter
export convergence_time
export mass_conservation
export check_iterConvergence
export check_timeConvergence
export check_timeConvergenceRelavite

include("CSIterConvergence.jl")

include("CSTimeConvergence.jl")

include("CSMassConservation.jl")

include("CSCheckConvergence.jl")

end # module
