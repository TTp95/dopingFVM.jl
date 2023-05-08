"""
dopingFVM

The module is structured in the following sub-modules:
- [`dopingFVM.Structures`](@ref)
- [`dopingFVM.StructuresConstructors`](@ref)
- [`dopingFVM.MeshGenerator`](@ref)
- [`dopingFVM.Tools`](@ref)
- [`dopingFVM.Bounds`](@ref)
- [`dopingFVM.Convergence`](@ref)
- [`dopingFVM.Diffusion`](@ref)
- [`dopingFVM.Convection`](@ref)
- [`dopingFVM.Source`](@ref)
- [`dopingFVM.Gradients`](@ref)
- [`dopingFVM.Transient`](@ref)
- [`dopingFVM.VelocityInterpolation`](@ref)
- [`dopingFVM.CoupleAlgorithms`](@ref)
- [`dopingFVM.Turbulence`](@ref)
- [`dopingFVM.TurbulenceCoupleAlgorithms`](@ref)
- [`dopingFVM.VisualizationData`](@ref)
- [`dopingFVM.Print`](@ref)
- [`dopingFVM.Solvers`](@ref)
- [`dopingFVM.dopingSolver`](@ref)
- [`dopingFVM.Backup`](@ref)
- [`dopingFVM.ScriptGenerator`](@ref)

The exported names are:
$(EXPORTS)
"""
module dopingFVM

using Base.Threads

using DocStringExtensions
using SparseArrays

include("Structures/Structures.jl")

include("StructuresConstructors/StructuresConstructors.jl")

include("MeshGenerator/MeshGenerator.jl")

include("Tools/Tools.jl")

include("Bounds/Bounds.jl")

include("Convergence/Convergence.jl")

include("Diffusion/Diffusion.jl")

include("Convection/Convection.jl")

include("Source/Source.jl")

include("Gradients/Gradients.jl")

include("Transient/Transient.jl")

include("VelocityInterpolation/VelocityInterpolation.jl")

include("CoupleAlgorithms/CoupleAlgorithms.jl")

include("Turbulence/Turbulence.jl")

include("TurbulenceCoupleAlgorithms/TurbulenceCoupleAlgorithms.jl")

include("VisualizationData/VisualizationData.jl")

include("Print/Print.jl")

include("Solvers/Solvers.jl")

include("dopingSolver/dopingSolver.jl")

include("Backup/Backup.jl")

include("ScriptGenerator/ScriptGenerator.jl")

include("Exports.jl")

end
