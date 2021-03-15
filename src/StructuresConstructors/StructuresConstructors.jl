"""
This module provides:
    -

The exported names are:
$(EXPORTS)
"""
module StructuresConstructors

using DocStringExtensions

using dopingFVM.Structures

export create_BoundStructured
export create_CSMaterial
export create_CSPhi
export create_CSFaceVelocity
export create_CSVelocity

include("BoundaryConditions.jl")

include("PhiCartesian.jl")

include("MaterialCartesian.jl")

include("PhiPolar.jl")

include("MaterialPolar.jl")

include("PhiCylindrical.jl")

include("MaterialCylindrical.jl")

include("PhiSpherical.jl")

include("MaterialSpherical.jl")

include("SystemConfig.jl")

include("SystemTime.jl")

include("DeltaTime.jl")

end # module
