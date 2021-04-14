"""

"""
function order_time! end

function order_time!(
    phi::UnionCSPhi,
)
    phi.time3 .= phi.time2
    phi.time2 .= phi.time1
    phi.time1 .= phi.eval

    return nothing
end

function order_time!(
    phi::CSVelocity1D,
)
    order_time!(phi.u)
    order_time!(phi.p)

    phi.fValues.uFaceTime3 .= phi.fValues.uFaceTime2
    phi.fValues.uFaceTime2 .= phi.fValues.uFaceTime1
    phi.fValues.uFaceTime1 .= phi.fValues.uFace

    return nothing
end

function order_time!(
    phi::CSVelocity2D,
)
    order_time!(phi.u)
    order_time!(phi.v)
    order_time!(phi.p)

    phi.fValues.uFaceTime3 .= phi.fValues.uFaceTime2
    phi.fValues.uFaceTime2 .= phi.fValues.uFaceTime1
    phi.fValues.uFaceTime1 .= phi.fValues.uFace

    phi.fValues.vFaceTime3 .= phi.fValues.vFaceTime2
    phi.fValues.vFaceTime2 .= phi.fValues.vFaceTime1
    phi.fValues.vFaceTime1 .= phi.fValues.vFace

    return nothing
end

function order_time!(
    phi::CSVelocity3D,
)
    order_time!(phi.u)
    order_time!(phi.v)
    order_time!(phi.w)
    order_time!(phi.p)

    phi.fValues.uFaceTime3 .= phi.fValues.uFaceTime2
    phi.fValues.uFaceTime2 .= phi.fValues.uFaceTime1
    phi.fValues.uFaceTime1 .= phi.fValues.uFace

    phi.fValues.vFaceTime3 .= phi.fValues.vFaceTime2
    phi.fValues.vFaceTime2 .= phi.fValues.vFaceTime1
    phi.fValues.vFaceTime1 .= phi.fValues.vFace

    phi.fValues.wFaceTime3 .= phi.fValues.wFaceTime2
    phi.fValues.wFaceTime2 .= phi.fValues.wFaceTime1
    phi.fValues.wFaceTime1 .= phi.fValues.wFace

    return nothing
end
