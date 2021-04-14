"""

"""
function order_iter! end

function order_iter!(
    phi::UnionCSPhi,
)
    phi.iter .= phi.eval

    return nothing
end

function order_iter!(
    phi::CSVelocity1D,
)
    order_iter!(phi.u)
    order_iter!(phi.p)

    phi.fValues.uFaceIter .= phi.fValues.uFace

    return nothing
end

function order_iter!(
    phi::CSVelocity2D,
)
    order_iter!(phi.u)
    order_iter!(phi.v)
    order_iter!(phi.p)

    phi.fValues.uFaceIter .= phi.fValues.uFace
    phi.fValues.vFaceIter .= phi.fValues.vFace

    return nothing
end

function order_iter!(
    phi::CSVelocity3D,
)
    order_iter!(phi.u)
    order_iter!(phi.v)
    order_iter!(phi.w)
    order_iter!(phi.p)

    phi.fValues.uFaceIter .= phi.fValues.uFace
    phi.fValues.vFaceIter .= phi.fValues.vFace
    phi.fValues.wFaceIter .= phi.fValues.wFace

    return nothing
end
