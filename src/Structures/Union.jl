"""


"""
UnionCSMesh1D = Union{CSMesh1D, CSMesh1DImmutable}

UnionCSMesh2D = Union{CSMesh2D, CSMesh2DImmutable}

UnionCSMesh3D = Union{CSMesh3D, CSMesh3DImmutable}

UnionCSMesh = Union{
    CSMesh1D,
    CSMesh1DImmutable,
    CSMesh2D,
    CSMesh2DImmutable,
    CSMesh3D,
    CSMesh3DImmutable
}

UnionCSPhi = Union{CSPhi1D, CSPhi2D, CSPhi3D}

UnionCSMaterial = Union{CSMaterial1D, CSMaterial2D, CSMaterial3D}

UnionCSConstantMaterial = Union{CSMaterialConstant, CSMaterialConstantImmutable}

UnionCSMaterialAll = Union{UnionCSMaterial, UnionCSConstantMaterial}
