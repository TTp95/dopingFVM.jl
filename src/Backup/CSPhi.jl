"""

"""
function save_Phi end

function save_Phi(
    phi::CSPhi1D,
    mesh::UnionCSMesh1D,
    name::String = "Phi1D",
)

    file_savePhi = open("./save_$(phi.key)_$(name).csv", "w")
    write(file_savePhi, "i,eval,iter,time1,time2,time3\n")

    for i in 1:mesh.l1
        write(file_savePhi, "$(i),$(phi.eval[i]),$(phi.iter[i]),$(phi.time1[i]),$(phi.time2[i]),$(phi.time3[i])\n")
    end

    file_savePhi = close(file_savePhi)

    return nothing
end

function save_Phi(
    phi::CSPhi2D,
    mesh::UnionCSMesh2D,
    name::String = "Phi2D",
)

    file_savePhi = open("./save_$(phi.key)_$(name).csv", "w")
    write(file_savePhi, "i,j,eval,iter,time1,time2,time3\n")

    for i in 1:mesh.l1
        for j in 1:mesh.m1
            write(file_savePhi, "$(i),$(j),$(phi.eval[i,j]),$(phi.iter[i,j]),$(phi.time1[i,j]),$(phi.time2[i,j]),$(phi.time3[i,j])\n")
        end
    end

    file_savePhi = close(file_savePhi)

    return nothing
end

function save_Phi(
    phi::CSPhi3D,
    mesh::UnionCSMesh3D,
    name::String = "Phi3D",
)

    file_savePhi = open("./save_$(phi.key)_$(name).csv", "w")
    write(file_savePhi, "i,j,k,eval,iter,time1,time2,time3\n")

    for i in 1:mesh.l1
        for j in 1:mesh.m1
            for k in 1:mesh.n1
                write(file_savePhi, "$(i),$(j),$(k),$(phi.eval[i,j,k]),$(phi.iter[i,j,k]),$(phi.time1[i,j,k]),$(phi.time2[i,j,k]),$(phi.time3[i,j,k])\n")
            end
        end
    end

    file_savePhi = close(file_savePhi)

    return nothing
end


"""

"""
function load_Phi! end

function load_Phi!(
    path::String,
    phi::CSPhi1D,
    mesh::UnionCSMesh1D,
)

    df = CSV.File(path; header=1, delim=",")

    if (df.i[end] != mesh.l1)
        error("")
    end

    last = length(df.i)

    for nn in 1:last
        phi.eval[df.i[nn]] = df.eval[nn]
        phi.iter[df.i[nn]] = df.iter[nn]
        phi.time1[df.i[nn]] = df.time1[nn]
        phi.time2[df.i[nn]] = df.time2[nn]
        phi.time3[df.i[nn]] = df.time3[nn]
    end


    return nothing
end

function load_Phi!(
    path::String,
    phi::CSPhi2D,
    mesh::UnionCSMesh2D,
)

    df = CSV.File(path; header=1, delim=",")

    if (df.i[end] != mesh.l1) || (df.j[end] != mesh.m1)
        error("")
    end

    last = length(df.i)

    for nn in 1:last
        phi.eval[df.i[nn],df.j[nn]] = df.eval[nn]
        phi.iter[df.i[nn],df.j[nn]] = df.iter[nn]
        phi.time1[df.i[nn],df.j[nn]] = df.time1[nn]
        phi.time2[df.i[nn],df.j[nn]] = df.time2[nn]
        phi.time3[df.i[nn],df.j[nn]] = df.time3[nn]
    end


    return nothing
end

function load_Phi!(
    path::String,
    phi::CSPhi3D,
    mesh::UnionCSMesh3D,
)

    df = CSV.File(path; header=1, delim=",")

    if (df.i[end] != mesh.l1) || (df.j[end] != mesh.m1) || (df.k[end] != mesh.n1)
        error("")
    end

    last = length(df.i)

    for nn in 1:last
        phi.eval[df.i[nn],df.j[nn],df.k[nn]] = df.eval[nn]
        phi.iter[df.i[nn],df.j[nn],df.k[nn]] = df.iter[nn]
        phi.time1[df.i[nn],df.j[nn],df.k[nn]] = df.time1[nn]
        phi.time2[df.i[nn],df.j[nn],df.k[nn]] = df.time2[nn]
        phi.time3[df.i[nn],df.j[nn],df.k[nn]] = df.time3[nn]
    end


    return nothing
end

