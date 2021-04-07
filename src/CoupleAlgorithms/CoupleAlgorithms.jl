"""
This module provides:
    -

The exported names are:
$(EXPORTS)
"""
module CoupleAlgorithms

using Base.Threads

using DocStringExtensions
using SparseArrays

using dopingFVM.Structures
using dopingFVM.Tools
using dopingFVM.Bounds
using dopingFVM.Gradients
using dopingFVM.Diffusion
using dopingFVM.Convection
using dopingFVM.Transient
using dopingFVM.VelocityInterpolation

export discretize_SIMPLE_PressureCorrection
export SIMPLE_correction!

export discretize_PPC_pressureEquation
export _projection_PC_pressure_laplacian_neighbors_
export _projection_PC_pressure_laplacian_
export _projection_PC_velocity_divergence_
export divergence_velocityToArray
export velocityProjection_PPC_Incremental!
export velocityProjection_PPC_nonIncremental!
export velocityProjection_PPC_Rotational!

include("./SIMPLE/CSPressureCorrectionDiscretization.jl")

include("./SIMPLE/CSCorrectionSIMPLE.jl")

include("./ProjectionPressureCorrection/CSPressureLaplacian.jl")

include("./ProjectionPressureCorrection/CSVelocityDivergence.jl")

include("./ProjectionPressureCorrection/CSVelocityDivergenceArray.jl")

include("./ProjectionPressureCorrection/CSPressureEquation.jl")

include("./ProjectionPressureCorrection/CSVelocityProjectionIncremental.jl")

include("./ProjectionPressureCorrection/CSVelocityProjectionNonIncremental.jl")

include("./ProjectionPressureCorrection/CSVelocityProjectionRotational.jl")

end # module
