"""
This module provides:
    -

The exported names are:
$(EXPORTS)
"""
module StructuresConstructors

using Base.Threads

using DocStringExtensions

using dopingFVM.Structures

export create_BoundsStructured
export create_Material
export create_Phi
export create_FaceVelocity
export create_Velocity
export crete_DeltaTime
export create_SystemConfig
export create_SystemControl
export create_SystemTime


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

include("SystemControl.jl")

include("SystemTime.jl")

include("DeltaTime.jl")

end # module
