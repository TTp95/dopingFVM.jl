"""

"""

function plot_paraviewVTK(
    system::SystemControl,
    mesh::MeshCartesianStructured,
    args...;
    fileVTK::Bool = false,
)
    n = length(args)

    if (typeof(mesh) <: UnionCSMesh1D)
        vtkfile = vtk_grid(
            ("./$(system.resultFolder)/paraview_$(system.problem)_$(system.timeSteps)"),
            mesh.x,
        )

    elseif (typeof(mesh) <: UnionCSMesh2D)
        vtkfile = vtk_grid(
            ("./$(system.resultFolder)/paraview_$(system.problem)_$(system.timeSteps)"),
            mesh.x,
            mesh.y,
        )

    elseif (typeof(mesh) <: UnionCSMesh3D)
        vtkfile = vtk_grid(
            ("./$(system.resultFolder)/paraview_$(system.problem)_$(system.timeSteps)"),
            mesh.x,
            mesh.y,
            mesh.z,
        )

    else
        error("Mesh type $(typeof(mesh)) incompatible...")
    end

    for x in 1:n
        if (typeof(args[x]) <: UnionCSVelocity)

            if (typeof(args[x]) <: CSVelocity1D)
                vtkfile["Velocity", VTKPointData()] = args[x].u.eval
                vtkfile["Pressure", VTKPointData()] = args[x].p.eval

            elseif (typeof(args[x]) <: CSVelocity2D)
                velocity = (args[x].u.eval, args[x].v.eval)
                vtkfile["Velocity", VTKPointData()] = velocity
                vtkfile["Pressure", VTKPointData()] = args[x].p.eval

            elseif (typeof(args[x]) <: CSVelocity3D)
                velocity = (args[x].u.eval, args[x].v.eval, args[x].w.eval)
                vtkfile["Velocity", VTKPointData()] = velocity
                vtkfile["Pressure", VTKPointData()] = args[x].p.eval
            end

        elseif (typeof(args[x]) <: PhiCartesianStructured)
            vtkfile["$(args[x].key)", VTKPointData()] =  args[x].eval

        elseif (typeof(args[x]) <: AbstractArray)
            vtkfile["array_$(x)", VTKPointData()] =  args[x]

        else
            error("Argumment number $(x + 2) incompatible...")
        end

    end

    if !fileVTK
        outfiles = vtk_save(vtkfile)
        return nothing
    else
        return vtkfile
    end
end
