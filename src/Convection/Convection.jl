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
using dopingFVM.Bounds

export discretize_convection

export _discretize_convection_centralDifference_
export _convection_centralDifference_neighbors_
export _convection_centralDifference_central_
export _convection_centralDifference_bounds_

export _discretize_convection_upwind_
export _convection_upwind_neighbors_
export _convection_upwind_central_
export _convection_upwind_bounds_

export _discretize_convection_downwind_
export _convection_downwind_neighbors_
export _convection_downwind_central_
export _convection_downwind_bounds_

export _discretize_convection_hybrid_
export _convection_hybrid_neighbors_
export _convection_hybrid_central_
export _convection_hybrid_bounds_

export _discretize_convection_powerlaw_
export _convection_powerlaw_neighbors_
export _convection_powerlaw_central_
export _convection_powerlaw_bounds_

export _discretize_convection_secondorderupwind_
export _convection_secondorderupwind_neighbors_
export _convection_secondorderupwind_central_
export _convection_secondorderupwind_bounds_

export _discretize_convection_quick_
export _convection_quick_neighbors_
export _convection_quick_central_
export _convection_quick_bounds_

export _discretize_convection_TVD_
export _discretize_convection_FROMM_

include("CSDiscretization.jl")

include("./CentralDifference/CSCentralDifference.jl")

include("./CentralDifference/CSCentralDifferenceCoefficients.jl")

include("./CentralDifference/CSCentralDifferenceCoefficientsBoundaries.jl")

include("./Upwind/CSUpwind.jl")

include("./Upwind/CSUpwindCoefficients.jl")

include("./Upwind/CSUpwindCoefficientsBoundaries.jl")

include("./Downwind/CSDownwind.jl")

include("./Downwind/CSDownwindCoefficients.jl")

include("./Downwind/CSDownwindCoefficientsBoundaries.jl")

include("./Hybrid/CSHybrid.jl")

include("./Hybrid/CSHybridCoefficients.jl")

include("./Hybrid/CSHybridCoefficientsBoundaries.jl")

include("./PowerLaw/CSPowerLaw.jl")

include("./PowerLaw/CSPowerLawCoefficients.jl")

include("./PowerLaw/CSPowerLawCoefficientsBoundaries.jl")

include("./SOU/CSSOU.jl")

include("./SOU/CSSOUCoefficients.jl")

include("./SOU/CSSOUCoefficientsBoundaries.jl")

include("./QUICK/CSQUICK.jl")

include("./QUICK/CSQUICKCoefficients.jl")

include("./QUICK/CSQUICKCoefficientsBoundaries.jl")

include("./TVD/TVD.jl")

include("./FROMM/FROMM.jl")

end # module
