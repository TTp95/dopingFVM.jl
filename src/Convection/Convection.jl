"""
This module provides:
    -

The exported names are:
$(EXPORTS)
"""
module Convection

using Base.Threads

using DocStringExtensions
using SparseArrays

using dopingFVM.Structures
using dopingFVM.Tools

export discretize_convection
export _discretize_convection_centralDifference_
export _convection_centralDifference_neighbors_
export _convection_centralDifference_central_
export _convection_centralDifference_bounds_

include("CSDiscretization.jl")

include("./CentralDifference/CSCentralDifference.jl")

include("./CentralDifference/CSCentralDifferenceCoefficients.jl")

include("./CentralDifference/CSCentralDifferenceCoefficientsBoundaries.jl")

end # module
