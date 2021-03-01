"""
dopingFVM

, grid-based approximation of PDEs in the Julia programming language

This module provides rich set of tools for the numerical solution of PDE, mainly based
on finite element methods.

The module is structured in the following sub-modules:
- [`dopingFVM.Structures`](@ref)

The exported names are:
$(EXPORTS)
"""
module dopingFVM

using DocStringExtensions

include("Convection/Convection.jl")

include("CoupleAlgorithms/CoupleAlgorithms.jl")

include("Diffusion/Diffusion.jl")

include("dopingSolver/dopingSolver.jl")

include("MeshGenerator/MeshGenerator.jl")

include("Source/Source.jl")

include("Structures/Structures.jl")

include("StructuresConstructors/StructuresConstructors.jl")

include("Tools/Tools.jl")

include("Transient/Transient.jl")

include("Turbulence/Turbulence.jl")

include("TurbulenceCoupleAlgorithms/TurbulenceCoupleAlgorithms.jl")

include("VisualizationData/VisualizationData.jl")

include("Exports.jl")

end
