"""
This module provides:
    -

The exported names are:
$(EXPORTS)
"""
module Transient

using Base.Threads

using DocStringExtensions
using SparseArrays

using dopingFVM.Structures
using dopingFVM.Tools
using dopingFVM.Bounds

export discretize_time
export discretize_euler_time
export discretize_crankNicolson_time
export discretize_BDF2_time

include("CSTransientTerm.jl")

include("CSEuler.jl")

include("CSCrankNicolson.jl")

include("CSBDF2.jl")

include("CSBDF3.jl")

end # module
