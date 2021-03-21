"""
This module provides:
    -

The exported names are:
$(EXPORTS)
"""
module Transient

using Base.Threads

using DocStringExtensions

using dopingFVM.Structures
using dopingFVM.Tools

export discretize_time

include("CSTransientTerm.jl")

end # module
