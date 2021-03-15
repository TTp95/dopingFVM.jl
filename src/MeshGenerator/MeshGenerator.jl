"""
This module provides:
    -

The exported names are:
$(EXPORTS)
"""
module MeshGenerator

using Base.Threads

using DocStringExtensions

using dopingFVM.Structures

export create_uniform_CSMesh
export create_nonuniform_CSMesh

export create_uniform_CSMesh1D
export create_uniform_CSMesh2D
export create_uniform_CSMesh3D

export create_nonuniform_CSMesh1D
export create_nonuniform_CSMesh2D
export create_nonuniform_CSMesh3D

include("MeshCartesianUniform.jl")

include("MeshCartesianNonUniform.jl")

include("MeshPolarUniform.jl")

include("MeshPolarNonUniform.jl")

include("MeshCylindricalUniform.jl")

include("MeshCylindricalNonUniform.jl")

include("MeshSphericalUniform.jl")

include("MeshSphericalNonUniform.jl")

end # module
