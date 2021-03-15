"""


"""
UnionCSMesh1D = Union{CSMesh1D, CSMesh1DImmutable}

UnionCSMesh2D = Union{CSMesh2D, CSMesh2DImmutable}

UnionCSMesh3D = Union{CSMesh3D, CSMesh3DImmutable}

UnionCSPhi = Union{CSPhi1D, CSPhi2D, CSPhi3D}

UnionCSMesh = Union{
    CSMesh1D,
    CSMesh1DImmutable,
    CSMesh2D,
    CSMesh2DImmutable,
    CSMesh3D,
    CSMesh3DImmutable
}

UnionConstantMaterial = Union{CSMaterialConstant, CSMaterialConstantImmutable}
