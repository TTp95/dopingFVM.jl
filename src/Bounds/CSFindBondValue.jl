"""

"""
function find_bondValue end

function find_bondValue(
    i::Signed,
    id::Signed,
    phi::CSPhi1D,
    mesh::UnionCSMesh1D,
    bounds::Dict{String,BoundsStructured},
    bcord::Char,
)

    boundvalue = phi.eval[i]

    if (bcord == 'w') || (bcord == 'W')
        p1 = 'w'
        p2 = 'W'
    elseif (bcord == 'e') || (bcord == 'E')
        p1 = 'e'
        p2 = 'E'
    end

    if phi.bounds[i]
        nbounds = phi.nbounds[i]
        for nn in  1:nbounds
            if (bounds["$(id)g$(nn)"].cord == p1) || (bounds["$(id)g$(nn)"].cord == p2)
                boundvalue = bounds["$(id)g$(nn)"].eval
                break
            end
        end
    else
        boundvalue = phi.eval[i]
    end

    return boundvalue
end

function find_bondValue(
    i::Signed,
    j::Signed,
    id::Signed,
    phi::CSPhi2D,
    mesh::UnionCSMesh2D,
    bounds::Dict{String,BoundsStructured},
    bcord::Char,
)

    boundvalue = phi.eval[i,j]

    if (bcord == 'w') || (bcord == 'W')
        p1 = 'w'
        p2 = 'W'
    elseif (bcord == 'e') || (bcord == 'E')
        p1 = 'e'
        p2 = 'E'
    elseif (bcord == 's') || (bcord == 'S')
        p1 = 's'
        p2 = 'S'
    elseif (bcord == 'n') || (bcord == 'N')
        p1 = 'n'
        p2 = 'N'
    end

    if phi.bounds[i,j]
        nbounds = phi.nbounds[i,j]
        for nn in  1:nbounds
            if (bounds["$(id)g$(nn)"].cord == p1) || (bounds["$(id)g$(nn)"].cord == p2)
                boundvalue = bounds["$(id)g$(nn)"].eval
                break
            end
        end
    else
        boundvalue = phi.eval[i,j]
    end

    return boundvalue
end

function find_bondValue(
    i::Signed,
    j::Signed,
    k::Signed,
    id::Signed,
    phi::CSPhi3D,
    mesh::UnionCSMesh3D,
    bounds::Dict{String,BoundsStructured},
    bcord::Char,
)

    boundvalue = phi.eval[i,j,k]

    if (bcord == 'w') || (bcord == 'w')
        p1 = 'w'
        p2 = 'W'
    elseif (bcord == 'e') || (bcord == 'E')
        p1 = 'e'
        p2 = 'E'
    elseif (bcord == 's') || (bcord == 'S')
        p1 = 's'
        p2 = 'S'
    elseif (bcord == 'n') || (bcord == 'N')
        p1 = 'n'
        p2 = 'N'
    elseif (bcord == 'b') || (bcord == 'B')
        p1 = 'b'
        p2 = 'B'
    elseif (bcord == 't') || (bcord == 'T')
        p1 = 't'
        p2 = 'T'
    end

    if phi.bounds[i,j,k]
        nbounds = phi.nbounds[i,j,k]
        for nn in  1:nbounds
            if (bounds["$(id)g$(nn)"].cord == p1) || (bounds["$(id)g$(nn)"].cord == p2)
                boundvalue = bounds["$(id)g$(nn)"].eval
                break
            end
        end
    else
        boundvalue = phi.eval[i,j,k]
    end

    return boundvalue
end
