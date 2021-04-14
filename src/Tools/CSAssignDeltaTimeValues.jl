"""

"""
function assign_dt!(
    deltat::DeltaTime,
    args...,
)
    nt = length(args)

    if (nt > 3)
        nt = 3
    end

    if (!deltat.transient)
        deltat.dt1 = 1.0e20
        deltat.dt2 = 1.0e20
        deltat.dt3 = 1.0e20
    elseif (nt == 1)
        deltat.dt1 = args[1]
        deltat.dt2 = args[1]
        deltat.dt3 = args[1]
    elseif (nt == 2)
        deltat.dt1 = args[1]
        deltat.dt2 = args[2]
        deltat.dt3 = args[2]
    elseif (nt == 3)
        deltat.dt1 = args[1]
        deltat.dt2 = args[2]
        deltat.dt3 = args[3]
    end

    return nothing
end

"""

"""
function push_dt!(
    deltat::DeltaTime,
    new::AbstractFloat,
)
    deltat.dt3 = deltat.dt2
    deltat.dt2 = deltat.dt1
    deltat.dt1 = new

    return nothing
end
