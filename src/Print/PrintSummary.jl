"""

"""
function print_summary(
    system::SystemControl,
    deltaTime::DeltaTime,
    timer::SystemTime,
    mesh::UnionCSMesh,
    #args..,
)
    summary = "$(system.resultFolder)/summary_$(system.problem).txt"
    file_summary = open(summary, "w")

    write(file_summary, "summary: $(system.problem) \n \n")

    print("\n")
    print("summary: $(system.problem) \n \n")

    if (typeof(mesh) <: UnionCSMesh1D)
        write(file_summary, "Mesh: $(mesh.l1) \n")
        write(file_summary, "Total de volumenes: $(mesh.l1) \n")
        print("Mesh: $(mesh.l1) \n")
        print("Total de volumenes: $(mesh.l1) \n")
    elseif (typeof(mesh) <: UnionCSMesh2D)
        write(file_summary, "Mesh: $(mesh.l1)x$(mesh.m1)  \n")
        write(file_summary, "Total de volumenes: $(mesh.l1 * mesh.m1) \n")
        print("Mesh: $(mesh.l1)x$(mesh.m1)  \n")
        print("Total de volumenes: $(mesh.l1 * mesh.m1) \n")
    elseif (typeof(mesh) <: UnionCSMesh3D)
        write(file_summary, "Mesh: $(mesh.l1)x$(mesh.m1)x$(mesh.n1) \n")
        write(file_summary, "Total de volumenes: $(mesh.l1 * mesh.m1 * mesh.n1) \n")
        print("Mesh: $(mesh.l1)x$(mesh.m1)x$(mesh.n1) \n")
        print("Total de volumenes: $(mesh.l1 * mesh.m1 * mesh.n1) \n")
    else
        error("Mesh type unimplemented")
    end

    write(file_summary, "\n")
    print("\n")

    if deltaTime.transient
        write(file_summary, "unsteady problem solution \n")
        write(file_summary, "\n")
        write(file_summary, "initial time: $(deltaTime.initialTime) \n")
        write(file_summary, "final time: $(deltaTime.finalTime) \n")
        write(file_summary, "uniform time step: $(deltaTime.dt1) \n")
        write(file_summary, "total time steps: $(system.timeSteps) \n \n")
        write(file_summary, "\n")
        print("unsteady problem solution \n")
        print("\n")
        print("initial time: $(deltaTime.initialTime) \n")
        print("final time: $(deltaTime.finalTime) \n")
        print("uniform time step: $(deltaTime.dt1) \n")
        print("total time steps: $(system.timeSteps) \n \n")
        print("\n")
    else
        write(file_summary, "steady state problem \n")
        write(file_summary, "\n")
        print("steady state problem \n")
        print("\n")
    end

    write(file_summary, "time CPU: $(timer.cpu) \n")
    write(file_summary, "total iterations: $(system.iterations) \n \n")
    write(file_summary, "time properties definition: $(timer.properties) \n")
    write(file_summary, "time discretization: $(timer.discretize) \n")
    write(file_summary, "time solvers: $(timer.solver) \n")
    write(file_summary, "time convergence check: $(timer.convergence) \n")
    write(file_summary, "time print results: $(timer.print) \n")
    write(file_summary, "\n")

    print("time CPU: $(timer.cpu) \n")
    print("total iterations: $(system.iterations) \n \n")
    print("time properties definition: $(timer.properties) \n")
    print("time discretization: $(timer.discretize) \n")
    print("time solvers: $(timer.solver) \n")
    print("time convergence check: $(timer.convergence) \n")
    print("time print results: $(timer.print) \n")
    print("\n")


    file_summary = close(file_summary)

    return nothing
end
