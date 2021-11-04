"""

"""
function save_System(
    system::SystemControl,
    dt::DeltaTime,
    name::String = "System",
)

    file_save = open("./save_$(name).csv", "w")
    write(file_save, "iterations,time,timeSteps,dt1,dt2,dt3\n")

    write(file_save, "$(system.iterations),$(system.time),$(system.timeSteps),$(dt.dt1),$(dt.dt2),$(dt.dt3)\n")

    file_save = close(file_save)

    return nothing
end

function load_System!(
    path::String,
    system::SystemControl,
    dt::DeltaTime,
)

    df = CSV.File(path; header=1, delim=",")

    system.iterations = df.iterations[end]
    system.time = df.time[end]
    system.timeSteps = df.timeSteps[end]

    dt.dt1 = df.dt1[end]
    dt.dt2 = df.dt2[end]
    dt.dt3 = df.dt3[end]

    return nothing
end

