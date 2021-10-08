"""

"""
function bounds_template end

function bounds_template(
    phi::CSPhi1D,
    mesh::UnionCSMesh1D;
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
)
    # vertices
    phi.bounds[1] = true
    phi.bounds[mesh.l1] = true
    phi.nbounds[1] = 1
    phi.nbounds[mesh.l1] = 1

    dict = create_BoundsDict(
        phi,
        mesh;
        T = T,
        N = N,
    )

    dict["$(phi.gIndex[1])g1"].cord = 'w'
    dict["$(phi.gIndex[mesh.l1])g1"].cord = 'e'

    return dict
end

function bounds_template(
    phi::CSPhi2D,
    mesh::UnionCSMesh2D;
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
)
    # horizontal aristas
    for i in 2:(mesh.l1 - 1)
        phi.bounds[i,1] = true
        phi.bounds[i,mesh.m1] = true
        phi.nbounds[i,1] = 1
        phi.nbounds[i,mesh.m1] = 1
    end

    # vertical aristas
    for j in 2:(mesh.m1 - 1)
        phi.bounds[1,j] = true
        phi.bounds[mesh.l1,j] = true
        phi.nbounds[1,j] = 1
        phi.nbounds[mesh.l1,j] = 1
    end

    phi.bounds[1,1] = true
    phi.nbounds[1,1] = 2
    phi.bounds[mesh.l1,mesh.m1] = true
    phi.nbounds[mesh.l1,mesh.m1] = 2
    phi.bounds[1,mesh.m1] = true
    phi.nbounds[1,mesh.m1] = 2
    phi.bounds[mesh.l1,1] = true
    phi.nbounds[mesh.l1,1] = 2

    dict = create_BoundsDict(
        phi,
        mesh;
        T = T,
        N = N,
    )

    for i in 2:(mesh.l1 - 1)
        dict["$(phi.gIndex[i,1])g1"].cord = 's'
        dict["$(phi.gIndex[i,mesh.m1])g1"].cord = 'n'
    end

    for j in 2:(mesh.m1 - 1)
        dict["$(phi.gIndex[1,j])g1"].cord = 'w'
        dict["$(phi.gIndex[mesh.l1,j])g1"].cord = 'e'
    end

    dict["$(phi.gIndex[1,1])g1"].cord = 'w'
    dict["$(phi.gIndex[1,1])g2"].cord = 's'
    dict["$(phi.gIndex[mesh.l1,mesh.m1])g1"].cord = 'e'
    dict["$(phi.gIndex[mesh.l1,mesh.m1])g2"].cord = 'n'
    dict["$(phi.gIndex[1,mesh.m1])g1"].cord = 'w'
    dict["$(phi.gIndex[1,mesh.m1])g2"].cord = 'n'
    dict["$(phi.gIndex[mesh.l1,1])g1"].cord = 'e'
    dict["$(phi.gIndex[mesh.l1,1])g2"].cord = 's'

    return dict
end

function bounds_template(
    phi::CSPhi3D,
    mesh::UnionCSMesh3D;
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
)
    # south and north surfaces
    for i in 2:(mesh.l1 - 1)
        for k in 2:(mesh.n1 - 1)
            phi.bounds[i,1,k] = true
            phi.bounds[i,mesh.m1,k] = true
            phi.nbounds[i,1,k] = 1
            phi.nbounds[i,mesh.m1,k] = 1
        end
    end

    # west and east surfaces
    for j in 2:(mesh.m1 - 1)
        for k in 2:(mesh.n1 - 1)
            phi.bounds[1,j,k] = true
            phi.bounds[mesh.l1,j,k] = true
            phi.nbounds[1,j,k] = 1
            phi.nbounds[mesh.l1,j,k] = 1
        end
    end

    # bottom and top surfaces
    for i in 2:(mesh.l1 - 1)
        for j in 2:(mesh.m1 - 1)
            phi.bounds[i,j,1] = true
            phi.bounds[i,j,mesh.n1] = true
            phi.nbounds[i,j,1] = 1
            phi.nbounds[i,j,mesh.n1] = 1
        end
    end

    # horizontal aristas
    for i in 2:(mesh.l1 - 1)
        phi.bounds[i,1,1] = true
        phi.bounds[i,mesh.m1,1] = true
        phi.bounds[i,1,mesh.n1] = true
        phi.bounds[i,mesh.m1,mesh.n1] = true
        phi.nbounds[i,1,1] = 2
        phi.nbounds[i,mesh.m1,1] = 2
        phi.nbounds[i,1,mesh.n1] = 2
        phi.nbounds[i,mesh.m1,mesh.n1] = 2
    end

    # vertical aristas
    for j in 2:(mesh.m1 - 1)
        phi.bounds[1,j,1] = true
        phi.bounds[mesh.l1,j,1] = true
        phi.bounds[1,j,mesh.n1] = true
        phi.bounds[mesh.l1,j,mesh.n1] = true
        phi.nbounds[1,j,1] = 2
        phi.nbounds[mesh.l1,j,1] = 2
        phi.nbounds[1,j,mesh.n1] = 2
        phi.nbounds[mesh.l1,j,mesh.n1] = 2
    end

    # deep aristas
    for k in 2:(mesh.n1 - 1)
        phi.bounds[1,1,k] = true
        phi.bounds[mesh.l1,1,k] = true
        phi.bounds[1,mesh.m1,k] = true
        phi.bounds[mesh.l1,mesh.m1,k] = true
        phi.nbounds[1,1,k] = 2
        phi.nbounds[mesh.l1,1,k] = 2
        phi.nbounds[1,mesh.m1,k] = 2
        phi.nbounds[mesh.l1,mesh.m1,k] = 2
    end

    # vertices
    phi.bounds[1,1,1] = true
    phi.bounds[mesh.l1,1,1] = true
    phi.bounds[mesh.l1,mesh.m1,1] = true
    phi.bounds[mesh.l1,1,mesh.n1] = true
    phi.bounds[1,mesh.m1,1] = true
    phi.bounds[1,mesh.m1,mesh.n1] = true
    phi.bounds[1,1,mesh.n1] = true
    phi.bounds[mesh.l1,mesh.m1,mesh.n1] = true
    phi.nbounds[1,1,1] = 3
    phi.nbounds[mesh.l1,1,1] = 3
    phi.nbounds[mesh.l1,mesh.m1,1] = 3
    phi.nbounds[mesh.l1,1,mesh.n1] = 3
    phi.nbounds[1,mesh.m1,1] = 3
    phi.nbounds[1,mesh.m1,mesh.n1] = 3
    phi.nbounds[1,1,mesh.n1] = 3
    phi.nbounds[mesh.l1,mesh.m1,mesh.n1] = 3

    dict = create_BoundsDict(
        phi,
        mesh;
        T = T,
        N = N,
    )

    # south and north surfaces
    for i in 2:(mesh.l1 - 1)
        for k in 2:(mesh.n1 - 1)
            dict["$(phi.gIndex[i,1,k])g1"].cord = 's'
            dict["$(phi.gIndex[i,mesh.m1,k])g1"].cord = 'n'
        end
    end

    # west and east surfaces
    for j in 2:(mesh.m1 - 1)
        for k in 2:(mesh.n1 - 1)
            dict["$(phi.gIndex[1,j,k])g1"].cord = 'w'
            dict["$(phi.gIndex[mesh.l1,j,k])g1"].cord = 'e'
        end
    end

    # bottom and top surfaces
    for i in 2:(mesh.l1 - 1)
        for j in 2:(mesh.m1 - 1)
            dict["$(phi.gIndex[i,j,1])g1"].cord = 'b'
            dict["$(phi.gIndex[i,j,mesh.n1])g1"].cord = 't'
        end
    end

    # horizontal aristas
    for i in 2:(mesh.l1 - 1)
        dict["$(phi.gIndex[i,1,1])g1"].cord = 's'
        dict["$(phi.gIndex[i,1,1])g2"].cord = 'b'
        dict["$(phi.gIndex[i,mesh.m1,1])g1"].cord = 'n'
        dict["$(phi.gIndex[i,mesh.m1,1])g2"].cord = 'b'
        dict["$(phi.gIndex[i,1,mesh.n1])g1"].cord = 's'
        dict["$(phi.gIndex[i,1,mesh.n1])g2"].cord = 't'
        dict["$(phi.gIndex[i,mesh.m1,mesh.n1])g1"].cord = 'n'
        dict["$(phi.gIndex[i,mesh.m1,mesh.n1])g2"].cord = 't'
    end

    # vertical aristas
    for j in 2:(mesh.m1 - 1)
        dict["$(phi.gIndex[1,j,1])g1"].cord = 'w'
        dict["$(phi.gIndex[1,j,1])g2"].cord = 'b'
        dict["$(phi.gIndex[mesh.l1,j,1])g1"].cord = 'e'
        dict["$(phi.gIndex[mesh.l1,j,1])g2"].cord = 'b'
        dict["$(phi.gIndex[1,j,mesh.n1])g1"].cord = 'w'
        dict["$(phi.gIndex[1,j,mesh.n1])g2"].cord = 't'
        dict["$(phi.gIndex[mesh.l1,j,mesh.n1])g1"].cord = 'e'
        dict["$(phi.gIndex[mesh.l1,j,mesh.n1])g2"].cord = 't'
    end

    # deep aristas
    for k in 2:(mesh.n1 - 1)
        dict["$(phi.gIndex[1,1,k])g1"].cord = 'w'
        dict["$(phi.gIndex[1,1,k])g2"].cord = 's'
        dict["$(phi.gIndex[mesh.l1,1,k])g1"].cord = 'e'
        dict["$(phi.gIndex[mesh.l1,1,k])g2"].cord = 's'
        dict["$(phi.gIndex[1,mesh.m1,k])g1"].cord = 'w'
        dict["$(phi.gIndex[1,mesh.m1,k])g2"].cord = 'n'
        dict["$(phi.gIndex[mesh.l1,mesh.m1,k])g1"].cord = 'e'
        dict["$(phi.gIndex[mesh.l1,mesh.m1,k])g2"].cord = 'n'
    end

    # vertices
    dict["$(phi.gIndex[1,1,1])g1"].cord = 'w'
    dict["$(phi.gIndex[1,1,1])g2"].cord = 's'
    dict["$(phi.gIndex[1,1,1])g3"].cord = 'b'
    dict["$(phi.gIndex[mesh.l1,1,1])g1"].cord = 'e'
    dict["$(phi.gIndex[mesh.l1,1,1])g2"].cord = 's'
    dict["$(phi.gIndex[mesh.l1,1,1])g3"].cord = 'b'
    dict["$(phi.gIndex[mesh.l1,mesh.m1,1])g1"].cord = 'e'
    dict["$(phi.gIndex[mesh.l1,mesh.m1,1])g2"].cord = 'n'
    dict["$(phi.gIndex[mesh.l1,mesh.m1,1])g3"].cord = 'b'
    dict["$(phi.gIndex[mesh.l1,1,mesh.n1])g1"].cord = 'e'
    dict["$(phi.gIndex[mesh.l1,1,mesh.n1])g2"].cord = 's'
    dict["$(phi.gIndex[mesh.l1,1,mesh.n1])g3"].cord = 't'
    dict["$(phi.gIndex[1,mesh.m1,1])g1"].cord = 'w'
    dict["$(phi.gIndex[1,mesh.m1,1])g2"].cord = 'n'
    dict["$(phi.gIndex[1,mesh.m1,1])g3"].cord = 'b'
    dict["$(phi.gIndex[1,mesh.m1,mesh.n1])g1"].cord = 'w'
    dict["$(phi.gIndex[1,mesh.m1,mesh.n1])g2"].cord = 'n'
    dict["$(phi.gIndex[1,mesh.m1,mesh.n1])g3"].cord = 't'
    dict["$(phi.gIndex[1,1,mesh.n1])g1"].cord = 'w'
    dict["$(phi.gIndex[1,1,mesh.n1])g2"].cord = 's'
    dict["$(phi.gIndex[1,1,mesh.n1])g3"].cord = 't'
    dict["$(phi.gIndex[mesh.l1,mesh.m1,mesh.n1])g1"].cord = 'e'
    dict["$(phi.gIndex[mesh.l1,mesh.m1,mesh.n1])g2"].cord = 'n'
    dict["$(phi.gIndex[mesh.l1,mesh.m1,mesh.n1])g3"].cord = 't'

    return dict
end
