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
