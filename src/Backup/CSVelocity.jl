"""

"""
function save_Velocity end

function save_Velocity(
    vel::CSVelocity1D,
    mesh::UnionCSMesh1D,
    name::String = "Velocity1D",
)

    # Velocity
    file_savePhi = open("./save_$(vel.u.key)_$(vel.p.key)_$(name).csv", "w")

    write(file_savePhi, "i,ueval,uiter,utime1,utime2,utime3,peval,piter,ptime1,ptime2,ptime3\n")

    for i in 1:mesh.l1
        write(file_savePhi, "$(i),$(vel.u.eval[i]),$(vel.u.iter[i]),$(vel.u.time1[i]),$(vel.u.time2[i]),$(vel.u.time3[i]),$(vel.p.eval[i]),$(vel.p.iter[i]),$(vel.p.time1[i]),$(vel.p.time2[i]),$(vel.p.time3[i])\n")
    end

    file_savePhi = close(file_savePhi)

    # Face Value U
    file_savePhi = open("./save_fValue_$(vel.u.key)_$(name).csv", "w")

    write(file_savePhi, "i,ueval,uiter,utime1,utime2,utime3\n")

    for i in 1:(mesh.l1+1)
        write(file_savePhi, "$(i),$(vel.fValues.uFace[i]),$(vel.fValues.uFaceIter[i]),$(vel.fValues.uFaceTime1[i]),$(vel.fValues.uFaceTime2[i,j]),$(vel.fValues.uFaceTime3[i])\n")
    end

    file_savePhi = close(file_savePhi)

    return nothing
end

function save_Velocity(
    vel::CSVelocity2D,
    mesh::UnionCSMesh2D,
    name::String = "Velocity2D",
)

    # Velocity
    file_savePhi = open("./save_$(vel.u.key)_$(vel.v.key)_$(vel.p.key)_$(name).csv", "w")

    write(file_savePhi, "i,j,ueval,uiter,utime1,utime2,utime3,veval,viter,vtime1,vtime2,vtime3,peval,piter,ptime1,ptime2,ptime3\n")

    for i in 1:mesh.l1
        for j in 1:mesh.m1
            write(file_savePhi, "$(i),$(j),$(vel.u.eval[i,j]),$(vel.u.iter[i,j]),$(vel.u.time1[i,j]),$(vel.u.time2[i,j]),$(vel.u.time3[i,j]),$(vel.v.eval[i,j]),$(vel.v.iter[i,j]),$(vel.v.time1[i,j]),$(vel.v.time2[i,j]),$(vel.v.time3[i,j]),$(vel.p.eval[i,j]),$(vel.p.iter[i,j]),$(vel.p.time1[i,j]),$(vel.p.time2[i,j]),$(vel.p.time3[i,j])\n")
        end
    end

    file_savePhi = close(file_savePhi)

    # Face Value U
    file_savePhi = open("./save_fValue_$(vel.u.key)_$(name).csv", "w")

    write(file_savePhi, "i,j,ueval,uiter,utime1,utime2,utime3\n")

    for i in 1:(mesh.l1+1)
        for j in 1:mesh.m1
            write(file_savePhi, "$(i),$(j),$(vel.fValues.uFace[i,j]),$(vel.fValues.uFaceIter[i,j]),$(vel.fValues.uFaceTime1[i,j]),$(vel.fValues.uFaceTime2[i,j]),$(vel.fValues.uFaceTime3[i,j])\n")
        end
    end

    file_savePhi = close(file_savePhi)

    # Face Value V
    file_savePhi = open("./save_fValue_$(vel.v.key)_$(name).csv", "w")

    write(file_savePhi, "i,j,veval,viter,vtime1,vtime2,vtime3\n")

    for i in 1:mesh.l1
        for j in (1:mesh.m1+1)
            write(file_savePhi, "$(i),$(j),$(vel.fValues.vFace[i,j]),$(vel.fValues.vFaceIter[i,j]),$(vel.fValues.vFaceTime1[i,j]),$(vel.fValues.vFaceTime2[i,j]),$(vel.fValues.vFaceTime3[i,j])\n")
        end
    end

    file_savePhi = close(file_savePhi)

    return nothing
end

function save_Velocity(
    vel::CSVelocity3D,
    mesh::UnionCSMesh3D,
    name::String = "Velocity3D",
)

    # Velocity
    file_savePhi = open("./save_$(vel.u.key)_$(vel.v.key)_$(vel.w.key)_$(vel.p.key)_$(name).csv", "w")

    write(file_savePhi, "i,j,k,ueval,uiter,utime1,utime2,utime3,veval,viter,vtime1,vtime2,vtime3,weval,witer,wtime1,wtime2,wtime3,peval,piter,ptime1,ptime2,ptime3\n")

    for i in 1:mesh.l1
        for j in 1:mesh.m1
            for k in 1:mesh.n1
                write(file_savePhi, "$(i),$(j),$(k),$(vel.u.eval[i,j,k]),$(vel.u.iter[i,j,k]),$(vel.u.time1[i,j,k]),$(vel.u.time2[i,j,k]),$(vel.u.time3[i,j,k]),$(vel.v.eval[i,j,k]),$(vel.v.iter[i,j,k]),$(vel.v.time1[i,j,k]),$(vel.v.time2[i,j,k]),$(vel.v.time3[i,j,k]),$(vel.w.eval[i,j,k]),$(vel.w.iter[i,j,k]),$(vel.w.time1[i,j,k]),$(vel.w.time2[i,j,k]),$(vel.w.time3[i,j,k]),$(vel.p.eval[i,j,k]),$(vel.p.iter[i,j,k]),$(vel.p.time1[i,j,k]),$(vel.p.time2[i,j,k]),$(vel.p.time3[i,j,k])\n")
            end
        end
    end

    file_savePhi = close(file_savePhi)

    # Face Value U
    file_savePhi = open("./save_fValue_$(vel.u.key)_$(name).csv", "w")

    write(file_savePhi, "i,j,k,ueval,uiter,utime1,utime2,utime3\n")

    for i in 1:(mesh.l1+1)
        for j in 1:mesh.m1
            for k in 1:mesh.n1
                write(file_savePhi, "$(i),$(j),$(k),$(vel.fValues.uFace[i,j,k]),$(vel.fValues.uFaceIter[i,j,k]),$(vel.fValues.uFaceTime1[i,j,k]),$(vel.fValues.uFaceTime2[i,j,k]),$(vel.fValues.uFaceTime3[i,j,k])\n")
            end
        end
    end

    file_savePhi = close(file_savePhi)

    # Face Value V
    file_savePhi = open("./save_fValue_$(vel.v.key)_$(name).csv", "w")

    write(file_savePhi, "i,j,k,veval,viter,vtime1,vtime2,vtime3\n")

    for i in 1:mesh.l1
        for j in 1:(mesh.m1+1)
            for k in 1:mesh.n1
                write(file_savePhi, "$(i),$(j),$(k),$(vel.fValues.vFace[i,j,k]),$(vel.fValues.vFaceIter[i,j,k]),$(vel.fValues.vFaceTime1[i,j,k]),$(vel.fValues.vFaceTime2[i,j,k]),$(vel.fValues.vFaceTime3[i,j,k])\n")
            end
        end
    end

    file_savePhi = close(file_savePhi)

    # Face Value W
    file_savePhi = open("./save_fValue_$(vel.w.key)_$(name).csv", "w")

    write(file_savePhi, "i,j,k,weval,witer,wtime1,wtime2,wtime3\n")

    for i in 1:mesh.l1
        for j in 1:mesh.m1
            for k in 1:(mesh.n1+1)
                write(file_savePhi, "$(i),$(j),$(k),$(vel.fValues.wFace[i,j,k]),$(vel.fValues.wFaceIter[i,j,k]),$(vel.fValues.wFaceTime1[i,j,k]),$(vel.fValues.wFaceTime2[i,j,k]),$(vel.fValues.wFaceTime3[i,j,k])\n")
            end
        end
    end

    file_savePhi = close(file_savePhi)

    return nothing
end

"""

"""
function load_Velocity! end

function load_Velocity!(
    vel::CSVelocity1D,
    mesh::UnionCSMesh1D,
    path1::String,
    path2::String,
)


    df1 = CSV.File(path1; header=1, delim=",")
    df2 = CSV.File(path2; header=1, delim=",")

    if (df1.i[end] != mesh.l1)
        error("")
    end

    last1 = length(df1.i)
    last2 = length(df2.i)

    for nn in 1:last1
        vel.u.eval[df1.i[nn]] = df1.ueval[nn]
        vel.u.iter[df1.i[nn]] = df1.uiter[nn]
        vel.u.time1[df1.i[nn]] = df1.utime1[nn]
        vel.u.time2[df1.i[nn]] = df1.utime2[nn]
        vel.u.time3[df1.i[nn]] = df1.utime3[nn]
        vel.p.eval[df1.i[nn]] = df1.peval[nn]
        vel.p.iter[df1.i[nn]] = df1.piter[nn]
        vel.p.time1[df1.i[nn]] = df1.ptime1[nn]
        vel.p.time2[df1.i[nn]] = df1.ptime2[nn]
        vel.p.time3[df1.i[nn]] = df1.ptime3[nn]
    end

    for nn in 1:last2
        vel.fValues.uFace[df2.i[nn]] = df2.ueval[nn]
        vel.fValues.uFaceIter[df2.i[nn]] = df2.uiter[nn]
        vel.fValues.uFaceTime1[df2.i[nn]] = df2.utime1[nn]
        vel.fValues.uFaceTime2[df2.i[nn]] = df2.utime2[nn]
        vel.fValues.uFaceTime3[df2.i[nn]] = df2.utime3[nn]
    end

    return nothing
end

function load_Velocity!(
    vel::CSVelocity2D,
    mesh::UnionCSMesh2D,
    path1::String,
    path2::String,
    path3::String,
)


    df1 = CSV.File(path1; header=1, delim=",")
    df2 = CSV.File(path2; header=1, delim=",")
    df3 = CSV.File(path3; header=1, delim=",")

    if (df1.i[end] != mesh.l1) || (df1.j[end] != mesh.m1)
        error("")
    end

    last1 = length(df1.i)
    last2 = length(df2.i)
    last3 = length(df3.i)

    for nn in 1:last1
        vel.u.eval[df1.i[nn],df1.j[nn]] = df1.ueval[nn]
        vel.u.iter[df1.i[nn],df1.j[nn]] = df1.uiter[nn]
        vel.u.time1[df1.i[nn],df1.j[nn]] = df1.utime1[nn]
        vel.u.time2[df1.i[nn],df1.j[nn]] = df1.utime2[nn]
        vel.u.time3[df1.i[nn],df1.j[nn]] = df1.utime3[nn]
        vel.v.eval[df1.i[nn],df1.j[nn]] = df1.veval[nn]
        vel.v.iter[df1.i[nn],df1.j[nn]] = df1.viter[nn]
        vel.v.time1[df1.i[nn],df1.j[nn]] = df1.vtime1[nn]
        vel.v.time2[df1.i[nn],df1.j[nn]] = df1.vtime2[nn]
        vel.v.time3[df1.i[nn],df1.j[nn]] = df1.vtime3[nn]
        vel.p.eval[df1.i[nn],df1.j[nn]] = df1.peval[nn]
        vel.p.iter[df1.i[nn],df1.j[nn]] = df1.piter[nn]
        vel.p.time1[df1.i[nn],df1.j[nn]] = df1.ptime1[nn]
        vel.p.time2[df1.i[nn],df1.j[nn]] = df1.ptime2[nn]
        vel.p.time3[df1.i[nn],df1.j[nn]] = df1.ptime3[nn]
    end

    for nn in 1:last2
        vel.fValues.uFace[df2.i[nn],df2.j[nn]] = df2.ueval[nn]
        vel.fValues.uFaceIter[df2.i[nn],df2.j[nn]] = df2.uiter[nn]
        vel.fValues.uFaceTime1[df2.i[nn],df2.j[nn]] = df2.utime1[nn]
        vel.fValues.uFaceTime2[df2.i[nn],df2.j[nn]] = df2.utime2[nn]
        vel.fValues.uFaceTime3[df2.i[nn],df2.j[nn]] = df2.utime3[nn]
    end

    for nn in 1:last3
        vel.fValues.vFace[df3.i[nn],df3.j[nn]] = df3.veval[nn]
        vel.fValues.vFaceIter[df3.i[nn],df3.j[nn]] = df3.viter[nn]
        vel.fValues.vFaceTime1[df3.i[nn],df3.j[nn]] = df3.vtime1[nn]
        vel.fValues.vFaceTime2[df3.i[nn],df3.j[nn]] = df3.vtime2[nn]
        vel.fValues.vFaceTime3[df3.i[nn],df3.j[nn]] = df3.vtime3[nn]
    end

    return nothing
end

function load_Velocity!(
    vel::CSVelocity3D,
    mesh::UnionCSMesh3D,
    path1::String,
    path2::String,
    path3::String,
    path4::String,
)


    df1 = CSV.File(path1; header=1, delim=",")
    df2 = CSV.File(path2; header=1, delim=",")
    df3 = CSV.File(path3; header=1, delim=",")
    df4 = CSV.File(path4; header=1, delim=",")

    if (df1.i[end] != mesh.l1) || (df1.j[end] != mesh.m1) || (df1.k[end] != mesh.n1)
        error("")
    end

    last1 = length(df1.i)
    last2 = length(df2.i)
    last3 = length(df3.i)
    last4 = length(df4.i)

    for nn in 1:last1
        vel.u.eval[df1.i[nn],df1.j[nn],df1.k[nn]] = df1.ueval[nn]
        vel.u.iter[df1.i[nn],df1.j[nn],df1.k[nn]] = df1.uiter[nn]
        vel.u.time1[df1.i[nn],df1.j[nn],df1.k[nn]] = df1.utime1[nn]
        vel.u.time2[df1.i[nn],df1.j[nn],df1.k[nn]] = df1.utime2[nn]
        vel.u.time3[df1.i[nn],df1.j[nn],df1.k[nn]] = df1.utime3[nn]
        vel.v.eval[df1.i[nn],df1.j[nn],df1.k[nn]] = df1.veval[nn]
        vel.v.iter[df1.i[nn],df1.j[nn],df1.k[nn]] = df1.viter[nn]
        vel.v.time1[df1.i[nn],df1.j[nn],df1.k[nn]] = df1.vtime1[nn]
        vel.v.time2[df1.i[nn],df1.j[nn],df1.k[nn]] = df1.vtime2[nn]
        vel.v.time3[df1.i[nn],df1.j[nn],df1.k[nn]] = df1.vtime3[nn]
        vel.w.eval[df1.i[nn],df1.j[nn],df1.k[nn]] = df1.weval[nn]
        vel.w.iter[df1.i[nn],df1.j[nn],df1.k[nn]] = df1.witer[nn]
        vel.w.time1[df1.i[nn],df1.j[nn],df1.k[nn]] = df1.wtime1[nn]
        vel.w.time2[df1.i[nn],df1.j[nn],df1.k[nn]] = df1.wtime2[nn]
        vel.w.time3[df1.i[nn],df1.j[nn],df1.k[nn]] = df1.wtime3[nn]
        vel.p.eval[df1.i[nn],df1.j[nn],df1.k[nn]] = df1.peval[nn]
        vel.p.iter[df1.i[nn],df1.j[nn],df1.k[nn]] = df1.piter[nn]
        vel.p.time1[df1.i[nn],df1.j[nn],df1.k[nn]] = df1.ptime1[nn]
        vel.p.time2[df1.i[nn],df1.j[nn],df1.k[nn]] = df1.ptime2[nn]
        vel.p.time3[df1.i[nn],df1.j[nn],df1.k[nn]] = df1.ptime3[nn]
    end

    for nn in 1:last2
        vel.fValues.uFace[df2.i[nn],df2.j[nn],df2.k[nn]] = df2.ueval[nn]
        vel.fValues.uFaceIter[df2.i[nn],df2.j[nn],df2.k[nn]] = df2.uiter[nn]
        vel.fValues.uFaceTime1[df2.i[nn],df2.j[nn],df2.k[nn]] = df2.utime1[nn]
        vel.fValues.uFaceTime2[df2.i[nn],df2.j[nn],df2.k[nn]] = df2.utime2[nn]
        vel.fValues.uFaceTime3[df2.i[nn],df2.j[nn],df2.k[nn]] = df2.utime3[nn]
    end

    for nn in 1:last3
        vel.fValues.vFace[df3.i[nn],df3.j[nn],df3.k[nn]] = df3.veval[nn]
        vel.fValues.vFaceIter[df3.i[nn],df3.j[nn],df3.k[nn]] = df3.viter[nn]
        vel.fValues.vFaceTime1[df3.i[nn],df3.j[nn],df3.k[nn]] = df3.vtime1[nn]
        vel.fValues.vFaceTime2[df3.i[nn],df3.j[nn],df3.k[nn]] = df3.vtime2[nn]
        vel.fValues.vFaceTime3[df3.i[nn],df3.j[nn],df3.k[nn]] = df3.vtime3[nn]
    end

    for nn in 1:last4
        vel.fValues.wFace[df4.i[nn],df4.j[nn],df4.k[nn]] = df4.weval[nn]
        vel.fValues.wFaceIter[df4.i[nn],df4.j[nn],df4.k[nn]] = df4.witer[nn]
        vel.fValues.wFaceTime1[df4.i[nn],df4.j[nn],df4.k[nn]] = df4.wtime1[nn]
        vel.fValues.wFaceTime2[df4.i[nn],df4.j[nn],df4.k[nn]] = df4.wtime2[nn]
        vel.fValues.wFaceTime3[df4.i[nn],df4.j[nn],df4.k[nn]] = df4.wtime3[nn]
    end

    return nothing
end
