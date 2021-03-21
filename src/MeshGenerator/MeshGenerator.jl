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

export create_uniform_Mesh
export create_nonuniform_Mesh

export create_uniform_Mesh1D
export create_uniform_Mesh2D
export create_uniform_Mesh3D

export create_nonuniform_Mesh1D
export create_nonuniform_Mesh2D
export create_nonuniform_Mesh3D

include("MeshCartesianUniform.jl")

include("MeshCartesianNonUniform.jl")

include("MeshPolarUniform.jl")

include("MeshPolarNonUniform.jl")

include("MeshCylindricalUniform.jl")

include("MeshCylindricalNonUniform.jl")

include("MeshSphericalUniform.jl")

include("MeshSphericalNonUniform.jl")

end # module
