"""
This module provides:
    -

The exported names are:
$(EXPORTS)
"""
module Structures

using DocStringExtensions

export TypeDopingFMV

export MeshStructured
export PhiStructured
export MaterialStructured
export BoundStructured
export DeltaTime

export MeshCartesianStructured
export PhiCartesianStructured
export MaterialCartesianStructured

export MeshPolarStructured
export PhiPolarStructured
export MaterialPolarStructured

export MeshCylindricalStructured
export PhiCylindricalStructured
export MaterialCylindricalStructured

export MeshSphericalStructured
export PhiSphericalStructured
export MaterialSphericalStructured

export CSMaterial1D
export CSMaterial2D
export CSMaterial3D
export CSMaterialConstant
export CSMaterialConstantImmutable

export CSMesh1D
export CSMesh2D
export CSMesh3D
export CSMesh1DImmutable
export CSMesh2DImmutable
export CSMesh3DImmutable

export CSPhi1D
export CSPhi2D
export CSPhi3D
export CSFaceVelocity1D
export CSFaceVelocity2D
export CSFaceVelocity3D
export CSVelocity1D
export CSVelocity2D
export CSVelocity3D

export UnionCSMesh1D
export UnionCSMesh2D
export UnionCSMesh3D
export UnionCSPhi
export UnionCSMesh
export UnionConstantMaterial


include("AbstractTypes.jl")

include("BoundaryConditions.jl")

include("DeltaTime.jl")

include("MeshCartesian.jl")

include("PhiCartesian.jl")

include("MaterialCartesian.jl")

include("MeshPolar.jl")

include("PhiPolar.jl")

include("MaterialPolar.jl")

include("MeshCylindrical.jl")

include("PhiCylindrical.jl")

include("MaterialCylindrical.jl")

include("MeshSpherical.jl")

include("PhiSpherical.jl")

include("MaterialSpherical.jl")

include("SystemConfig.jl")

include("SystemTime.jl")

include("Union.jl")

end # module
