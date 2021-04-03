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

    boundvalue = 0.0

    if (bcord == 'w')
        p1 = 'w'
        p2 = 'W'
    elseif (bcord == 'e')
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
        boundvalue = 0.0
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

    boundvalue = 0.0

    if (bcord == 'w')
        p1 = 'w'
        p2 = 'W'
    elseif (bcord == 'e')
        p1 = 'e'
        p2 = 'E'
    elseif (bcord == 's')
        p1 = 's'
        p2 = 'S'
    elseif (bcord == 'n')
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
        boundvalue = 0.0
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

    boundvalue = 0.0

    if (bcord == 'w')
        p1 = 'w'
        p2 = 'W'
    elseif (bcord == 'e')
        p1 = 'e'
        p2 = 'E'
    elseif (bcord == 's')
        p1 = 's'
        p2 = 'S'
    elseif (bcord == 'n')
        p1 = 'n'
        p2 = 'N'
    elseif (bcord == 'b')
        p1 = 'b'
        p2 = 'B'
    elseif (bcord == 't')
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
        boundvalue = 0.0
    end

    return boundvalue
end
