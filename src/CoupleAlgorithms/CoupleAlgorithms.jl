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
using dopingFVM.Gradients
using dopingFVM.Diffusion
using dopingFVM.Convection
using dopingFVM.Transient
using dopingFVM.VelocityInterpolation

export discretize_SIMPLE_PressureCorrection
export _SIMPLE_PressureCorrection_Coefficients_
export SIMPLE_correction!

include("./SIMPLE/CSPressureCorrectionDiscretization.jl")

include("./SIMPLE/CSPressureCorrectionCoefficients.jl")

include("./SIMPLE/CSCorrectionSIMPLE.jl")

end # module
