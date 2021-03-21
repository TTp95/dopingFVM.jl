"""
This module provides:
    -

The exported names are:
$(EXPORTS)
"""
module Diffusion

using Base.Threads

using DocStringExtensions
using SparseArrays

using dopingFVM.Structures
using dopingFVM.Tools

export discretize_diffusion
export _discretize_diffusion_centralDifference_
export _diffusion_centralDifference_neighbors_
export _diffusion_centralDifference_central_
export _diffusion_centralDifference_bounds_

include("CSDiscretization.jl")

include("./CentralDifference/CSCentralDifference.jl")

include("./CentralDifference/CSCentralDifferenceCoefficients.jl")

include("./CentralDifference/CSCentralDifferenceCoefficientsBoundaries.jl")

end # module
