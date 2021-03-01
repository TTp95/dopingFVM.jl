"""
TypeDopingFMV is an abstract type that encompasses all the 'dopingFVM' structs.
"""
abstract type TypeDopingFMV end

# ---------------------------------------------------- #
"""
MeshCartesianStructured is an abstract type for al the cartesian structured meshes.
"""
abstract type MeshStructured <: TypeDopingFMV end

"""
MeshCartesianStructured is an abstract type for al the cartesian structured meshes.
"""
abstract type PhiStructured <: TypeDopingFMV end

"""
MeshCartesianStructured is an abstract type for al the cartesian structured meshes.
"""
abstract type MaterialStructured <: TypeDopingFMV end

# ---------------------------------------------------- #
"""
MeshCartesianStructured is an abstract type for al the cartesian structured meshes.
"""
abstract type MeshCartesianStructured <: MeshStructured end

"""
PhiCartesianStructured is an abstract type for al the cartesian structured Phi.
"""
abstract type PhiCartesianStructured <: PhiStructured end

"""
PhiCartesianStructured is an abstract type for al the cartesian structured Phi.
"""
abstract type MaterialCartesianStructured <: MaterialStructured end

# ---------------------------------------------------- #
"""
MeshCartesianStructured is an abstract type for al the cartesian structured meshes.
"""
abstract type MeshPolarStructured <: MeshStructured end

"""
PhiCartesianStructured is an abstract type for al the cartesian structured Phi.
"""
abstract type PhiPolarStructured <: PhiStructured end

"""
PhiCartesianStructured is an abstract type for al the cartesian structured Phi.
"""
abstract type MaterialPolarStructured <: MaterialStructured end

# ---------------------------------------------------- #
"""
MeshCartesianStructured is an abstract type for al the cartesian structured meshes.
"""
abstract type MeshCylindricalStructured <: MeshStructured end

"""
PhiCartesianStructured is an abstract type for al the cartesian structured Phi.
"""
abstract type PhiCylindricalStructured <: PhiStructured end

"""
PhiCartesianStructured is an abstract type for al the cartesian structured Phi.
"""
abstract type MaterialCylindricalStructured <: MaterialStructured end

# ---------------------------------------------------- #
"""
MeshCartesianStructured is an abstract type for al the cartesian structured meshes.
"""
abstract type MeshSphericalStructured <: MeshStructured end

"""
PhiCartesianStructured is an abstract type for al the cartesian structured Phi.
"""
abstract type PhiSphericalStructured <: PhiStructured end

"""
PhiCartesianStructured is an abstract type for al the cartesian structured Phi.
"""
abstract type MaterialSphericalStructured <: MaterialStructured end
