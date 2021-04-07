"""

"""
function SystemTime_cpu(
    Stime::SystemTime,
    switch::Bool;
)
    if switch
        Stime.cpuPoint = time()
    else
        Stime.cpu += time() - Stime.cpuPoint
        Stime.cpuPoint = 0.0
    end

    return nothing
end

"""

"""
function SystemTime_properties(
    Stime::SystemTime,
    switch::Bool;
)
    if switch
        Stime.propertiesPoint = time()
    else
        Stime.properties += time() - Stime.propertiesPoint
        Stime.propertiesPoint = 0.0
    end

    return nothing
end

"""

"""
function SystemTime_discretize(
    Stime::SystemTime,
    switch::Bool;
)
    if switch
        Stime.discretizePoint = time()
    else
        Stime.discretize += time() - Stime.discretizePoint
        Stime.discretizePoint = 0.0
    end

    return nothing
end

"""

"""
function SystemTime_solver(
    Stime::SystemTime,
    switch::Bool;
)
    if switch
        Stime.solverPoint = time()
    else
        Stime.solver += time() - Stime.solverPoint
        Stime.solverPoint = 0.0
    end

    return nothing
end

"""

"""
function SystemTime_convergence(
    Stime::SystemTime,
    switch::Bool;
)
    if switch
        Stime.convergencePoint = time()
    else
        Stime.convergence += time() - Stime.convergencePoint
        Stime.convergencePoint = 0.0
    end

    return nothing
end

"""

"""
function SystemTime_print(
    Stime::SystemTime,
    switch::Bool;
)
    if switch
        Stime.printPoint = time()
    else
        Stime.print += time() - Stime.printPoint
        Stime.printPoint = 0.0
    end

    return nothing
end
