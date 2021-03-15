"""
dopingFVM

,grid-based approximation of PDEs in the Julia programming language

This module provides rich set of tools for the numerical solution of PDE, mainly based
on finite element methods.

The module is structured in the following sub-modules:
- [`dopingFVM.Structures`](@ref)
- [`dopingFVM.StructuresConstructors`](@ref)
- [`dopingFVM.Tools`](@ref)
- [`dopingFVM.MeshGenerator`](@ref)
- [`dopingFVM.Convergence`](@ref)
- [`dopingFVM.Diffusion`](@ref)
- [`dopingFVM.Convection`](@ref)
- [`dopingFVM.Source`](@ref)
- [`dopingFVM.Transient`](@ref)
- [`dopingFVM.CoupleAlgorithms`](@ref)
- [`dopingFVM.Turbulence`](@ref)
- [`dopingFVM.TurbulenceCoupleAlgorithms`](@ref)
- [`dopingFVM.VisualizationData`](@ref)
- [`dopingFVM.dopingSolver`](@ref)
- [`dopingFVM.ScriptGenerator`](@ref)

The exported names are:
$(EXPORTS)
"""
module dopingFVM

using Base.Threads

using DocStringExtensions

include("Structures/Structures.jl")

include("StructuresConstructors/StructuresConstructors.jl")

include("MeshGenerator/MeshGenerator.jl")

include("Tools/Tools.jl")

include("Convergence/Convergence.jl")

include("Diffusion/Diffusion.jl")

include("Convection/Convection.jl")

include("Source/Source.jl")

include("Transient/Transient.jl")

include("CoupleAlgorithms/CoupleAlgorithms.jl")

include("Turbulence/Turbulence.jl")

include("TurbulenceCoupleAlgorithms/TurbulenceCoupleAlgorithms.jl")

include("VisualizationData/VisualizationData.jl")

include("dopingSolver/dopingSolver.jl")

include("ScriptGenerator/ScriptGenerator.jl")

include("Exports.jl")

end
