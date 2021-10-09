"""

"""
function compute_RhieChow_Time end

function compute_RhieChow_Time(
    AU::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,2},
    },
    velocity::CSVelocity1D,
    mesh::UnionCSMesh1D,
    deltat::DeltaTime,
    material::CSMaterial1D;
    materialtime1::CSMaterial1D = material,
    materialtime2::CSMaterial1D = material,
    materialtime3::CSMaterial1D = material,
    T::Type{<:AbstractFloat} = Float64,
    scheme::Signed = 1,
)
    rhieChow_Time = zeros(T, mesh.l1)

    if mesh.l1 != 1
        for i in 2:mesh.l1
            if velocity.u.onoff[i-1] && velocity.u.onoff[i]
                dx1 = 0.0
                dx2 =  0.0
                id1 = 0
                id2 = 0
                D1  = 0.0
                D2  = 0.0
                Df = 0.0
                Vf = 0.0
                coef1 = 0.0
                coef2 = 0.0
                coef3 = 0.0
                a = 0.0
                b = 0.0
                meantime1 = 0.0
                meantime2 = 0.0
                meantime3 = 0.0
                fValue_old1 = 0.0
                fValue_old2 = 0.0
                fValue_old3 = 0.0
                term1 = 0.0
                term2 = 0.0
                term3 = 0.0

                #parameters
                dx1 = 0.5 * mesh.dx[i-1]
                dx2 =  0.5 * mesh.dx[i]
                id1 = velocity.u.gIndex[i-1]
                id2 = velocity.u.gIndex[i]

                D1  = mesh.vol[i-1] / AU[id1,id1]
                D2  = mesh.vol[i] / AU[id2,id2]
                num = D1 * dx2 + D2 * dx1
                den = dx1 + dx2
                Df = num / den

                num = (0.5 * mesh.dx[i-1]) * dx2 + (0.5 * mesh.dx[i]) * dx1
                den = dx1 + dx2
                Vf = num / den

                if (scheme == 1) #euler
                    num = materialtime1.ρ[i-1] * mesh.vol[i-1]
                    den = deltat.dt1
                    a = num / den

                    num = materialtime1.ρ[i] * mesh.vol[i]
                    den = deltat.dt1
                    b  = num / den

                    num = a * dx2 + b * dx1
                    den = dx1 + dx2
                    coef1 =  num / den

                    num = velocity.u.time1[i-1] * dx2 + velocity.u.time1[i] * dx1
                    den = dx1 + dx2
                    meantime1 = num / den

                    fValue_old1 = velocity.fValues.uFaceTime1[i]

                    rhieChow_Time[i] = (coef1 * Df / Vf) * (-1.0) * (fValue_old1 - meantime1)

                elseif (scheme == 2) #CN
                    a1 = ( deltat.dt2 / (deltat.dt1 * (deltat.dt1 + deltat.dt2)))
                    a2 = ( (deltat.dt2 - deltat.dt1) / (deltat.dt1 + deltat.dt2) )
                    a3 = ( deltat.dt1 / (deltat.dt2 * (deltat.dt1 + deltat.dt2)))

                    #coef1
                    a = a2 * materialtime1.ρ[i-1] * mesh.vol[i-1]
                    b  = a2 * materialtime1.ρ[i] * mesh.vol[i]

                    num = a * dx2 + b * dx1
                    den = dx1 + dx2
                    coef1 =  num / den

                    num = velocity.u.time1[i-1] * dx2 + velocity.u.time1[i] * dx1
                    den = dx1 + dx2
                    meantime1 = num / den

                    fValue_old1 = velocity.fValues.uFaceTime1[i]

                    term1 = (coef1 * Df / Vf) * (-1.0) * (fValue_old1 - meantime1)

                    #coef2
                    a = a3 * materialtime2.ρ[i-1] * mesh.vol[i-1]
                    b  = a3 * materialtime2.ρ[i] * mesh.vol[i]

                    num = a * dx2 + b * dx1
                    den = dx1 + dx2
                    coef2 =  num / den

                    num = velocity.u.time2[i-1] * dx2 + velocity.u.time2[i] * dx1
                    den = dx1 + dx2
                    meantime2 = num / den

                    fValue_old2 = velocity.fValues.uFaceTime2[i]

                    term2 = (coef2 * Df / Vf) * (-1.0) * (fValue_old2 - meantime2)

                    #Final
                    rhieChow_Time[i] = term1 + term2

                elseif (scheme == 3) #BDF2
                    a1 = ( 1.0 / deltat.dt1 ) + ( 1.0 / (deltat.dt1 + deltat.dt2))
                    a2 = (( 1.0 / deltat.dt1 ) + ( 1.0 / deltat.dt2))
                    a3 = -1.0 * ( deltat.dt1 / (deltat.dt2 * (deltat.dt1 + deltat.dt2)))

                    #coef1
                    a = a2 * materialtime1.ρ[i-1] * mesh.vol[i-1]
                    b  = a2 * materialtime1.ρ[i] * mesh.vol[i]

                    num = a * dx2 + b * dx1
                    den = dx1 + dx2
                    coef1 =  num / den

                    num = velocity.u.time1[i-1] * dx2 + velocity.u.time1[i] * dx1
                    den = dx1 + dx2
                    meantime1 = num / den

                    fValue_old1 = velocity.fValues.uFaceTime1[i]

                    term1 = (coef1 * Df / Vf) * (-1.0) * (fValue_old1 - meantime1)

                    #coef2
                    a = a3 * materialtime2.ρ[i-1] * mesh.vol[i-1]
                    b  = a3 * materialtime2.ρ[i] * mesh.vol[i]

                    num = a * dx2 + b * dx1
                    den = dx1 + dx2
                    coef2 =  num / den

                    num = velocity.u.time2[i-1] * dx2 + velocity.u.time2[i] * dx1
                    den = dx1 + dx2
                    meantime2 = num / den

                    fValue_old2 = velocity.fValues.uFaceTime2[i]

                    term2 = (coef2 * Df / Vf) * (-1.0) * (fValue_old2 - meantime2)

                    #Final
                    rhieChow_Time[i] = term1 + term2

                else
                    error("Time scheme: $(scheme) unimplemented")

                end

            end

        end

    end

    return rhieChow_Time
end

function compute_RhieChow_Time(
    AU::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,2},
    },
    AV::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,2},
    },
    velocity::CSVelocity2D,
    mesh::UnionCSMesh2D,
    deltat::DeltaTime,
    material::CSMaterial2D;
    materialtime1::CSMaterial2D = material,
    materialtime2::CSMaterial2D = material,
    materialtime3::CSMaterial2D = material,
    T::Type{<:AbstractFloat} = Float64,
    scheme::Signed = 1,
)
    u_rhieChow_Time = zeros(T, mesh.l1, mesh.m1)
    v_rhieChow_Time = zeros(T, mesh.l1, mesh.m1)

    #Compute u
    if mesh.l1 != 1
        for i in 2:mesh.l1
            for j in 1:mesh.m1
                if velocity.u.onoff[i-1,j] && velocity.u.onoff[i,j]
                    dx1 = 0.0
                    dx2 =  0.0
                    id1 = 0
                    id2 = 0
                    D1  = 0.0
                    D2  = 0.0
                    Df = 0.0
                    Vf = 0.0
                    coef1 = 0.0
                    coef2 = 0.0
                    coef3 = 0.0
                    a = 0.0
                    b = 0.0
                    meantime1 = 0.0
                    meantime2 = 0.0
                    meantime3 = 0.0
                    fValue_old1 = 0.0
                    fValue_old2 = 0.0
                    fValue_old3 = 0.0
                    term1 = 0.0
                    term2 = 0.0
                    term3 = 0.0

                    #parameters
                    dx1 = 0.5 * mesh.dx[i-1]
                    dx2 =  0.5 * mesh.dx[i]
                    id1 = velocity.u.gIndex[i-1,j]
                    id2 = velocity.u.gIndex[i,j]

                    D1  = mesh.vol[i-1,j] / AU[id1,id1]
                    D2  = mesh.vol[i,j] / AU[id2,id2]
                    num = D1 * dx2 + D2 * dx1
                    den = dx1 + dx2
                    Df = num / den

                    num = (0.5 * mesh.dx[i-1]) * dx2 + (0.5 * mesh.dx[i]) * dx1
                    den = dx1 + dx2
                    Vf = num / den

                    if (scheme == 1) #euler
                        num = materialtime1.ρ[i-1,j] * mesh.vol[i-1,j]
                        den = deltat.dt1
                        a = num / den

                        num = materialtime1.ρ[i,j] * mesh.vol[i,j]
                        den = deltat.dt1
                        b  = num / den

                        num = a * dx2 + b * dx1
                        den = dx1 + dx2
                        coef1 =  num / den

                        num = velocity.u.time1[i-1,j] * dx2 + velocity.u.time1[i,j] * dx1
                        den = dx1 + dx2
                        meantime1 = num / den

                        fValue_old1 = velocity.fValues.uFaceTime1[i,j]

                        u_rhieChow_Time[i,j] = (coef1 * Df / Vf) * (-1.0) * (fValue_old1 - meantime1)

                    elseif (scheme == 2) #CN
                        a1 = ( deltat.dt2 / (deltat.dt1 * (deltat.dt1 + deltat.dt2)))
                        a2 = ( (deltat.dt2 - deltat.dt1) / (deltat.dt1 + deltat.dt2) )
                        a3 = ( deltat.dt1 / (deltat.dt2 * (deltat.dt1 + deltat.dt2)))

                        #coef1
                        a = a2 * materialtime1.ρ[i-1,j] * mesh.vol[i-1,j]
                        b  = a2 * materialtime1.ρ[i,j] * mesh.vol[i,j]

                        num = a * dx2 + b * dx1
                        den = dx1 + dx2
                        coef1 =  num / den

                        num = velocity.u.time1[i-1,j] * dx2 + velocity.u.time1[i,j] * dx1
                        den = dx1 + dx2
                        meantime1 = num / den

                        fValue_old1 = velocity.fValues.uFaceTime1[i,j]

                        term1 = (coef1 * Df / Vf) * (-1.0) * (fValue_old1 - meantime1)

                        #coef2
                        a = a3 * materialtime2.ρ[i-1,j] * mesh.vol[i-1,j]
                        b  = a3 * materialtime2.ρ[i,j] * mesh.vol[i,j]

                        num = a * dx2 + b * dx1
                        den = dx1 + dx2
                        coef2 =  num / den

                        num = velocity.u.time2[i-1,j] * dx2 + velocity.u.time2[i,j] * dx1
                        den = dx1 + dx2
                        meantime2 = num / den

                        fValue_old2 = velocity.fValues.uFaceTime2[i,j]

                        term2 = (coef2 * Df / Vf) * (-1.0) * (fValue_old2 - meantime2)

                        #Final
                        u_rhieChow_Time[i,j] = term1 + term2

                    elseif (scheme == 3) #BDF2
                        a1 = ( 1.0 / deltat.dt1 ) + ( 1.0 / (deltat.dt1 + deltat.dt2))
                        a2 = (( 1.0 / deltat.dt1 ) + ( 1.0 / deltat.dt2))
                        a3 = -1.0 * ( deltat.dt1 / (deltat.dt2 * (deltat.dt1 + deltat.dt2)))

                        #coef1
                        a = a2 * materialtime1.ρ[i-1,j] * mesh.vol[i-1,j]
                        b  = a2 * materialtime1.ρ[i,j] * mesh.vol[i,j]

                        num = a * dx2 + b * dx1
                        den = dx1 + dx2
                        coef1 =  num / den

                        num = velocity.u.time1[i-1,j] * dx2 + velocity.u.time1[i,j] * dx1
                        den = dx1 + dx2
                        meantime1 = num / den

                        fValue_old1 = velocity.fValues.uFaceTime1[i,j]

                        term1 = (coef1 * Df / Vf) * (-1.0) * (fValue_old1 - meantime1)

                        #coef2
                        a = a3 * materialtime2.ρ[i-1,j] * mesh.vol[i-1,j]
                        b  = a3 * materialtime2.ρ[i,j] * mesh.vol[i,j]

                        num = a * dx2 + b * dx1
                        den = dx1 + dx2
                        coef2 =  num / den

                        num = velocity.u.time2[i-1,j] * dx2 + velocity.u.time2[i,j] * dx1
                        den = dx1 + dx2
                        meantime2 = num / den

                        fValue_old2 = velocity.fValues.uFaceTime2[i,j]

                        term2 = (coef2 * Df / Vf) * (-1.0) * (fValue_old2 - meantime2)

                        #Final
                        u_rhieChow_Time[i,j] = term1 + term2

                    else
                        error("Time scheme: $(scheme) unimplemented")

                    end

                end
            end
        end
    end

    #Compute v
    if mesh.m1 != 1
        for i in 1:mesh.l1
            for j in 2:mesh.m1
                if velocity.v.onoff[i,j-1] && velocity.v.onoff[i,j]
                    dx1 = 0.0
                    dx2 =  0.0
                    id1 = 0
                    id2 = 0
                    D1  = 0.0
                    D2  = 0.0
                    Df = 0.0
                    Vf = 0.0
                    coef1 = 0.0
                    coef2 = 0.0
                    coef3 = 0.0
                    a = 0.0
                    b = 0.0
                    meantime1 = 0.0
                    meantime2 = 0.0
                    meantime3 = 0.0
                    fValue_old1 = 0.0
                    fValue_old2 = 0.0
                    fValue_old3 = 0.0
                    term1 = 0.0
                    term2 = 0.0
                    term3 = 0.0

                    #parameters
                    dx1 = 0.5 * mesh.dy[j-1]
                    dx2 =  0.5 * mesh.dy[j]
                    id1 = velocity.v.gIndex[i,j-1]
                    id2 = velocity.v.gIndex[i,j]

                    D1  = mesh.vol[i,j-1] / AV[id1,id1]
                    D2  = mesh.vol[i,j] / AV[id2,id2]
                    num = D1 * dx2 + D2 * dx1
                    den = dx1 + dx2
                    Df = num / den

                    num = (0.5 * mesh.dy[j-1]) * dx2 + (0.5 * mesh.dy[j]) * dx1
                    den = dx1 + dx2
                    Vf = num / den

                    if (scheme == 1) #euler
                        num = materialtime1.ρ[i,j-1] * mesh.vol[i,j-1]
                        den = deltat.dt1
                        a = num / den

                        num = materialtime1.ρ[i,j] * mesh.vol[i,j]
                        den = deltat.dt1
                        b  = num / den

                        num = a * dx2 + b * dx1
                        den = dx1 + dx2
                        coef1 =  num / den

                        num = velocity.v.time1[i,j-1] * dx2 + velocity.v.time1[i,j] * dx1
                        den = dx1 + dx2
                        meantime1 = num / den

                        fValue_old1 = velocity.fValues.vFaceTime1[i,j]

                        v_rhieChow_Time[i,j] = (coef1 * Df / Vf) * (-1.0) * (fValue_old1 - meantime1)

                    elseif (scheme == 2) #CN
                        a1 = ( deltat.dt2 / (deltat.dt1 * (deltat.dt1 + deltat.dt2)))
                        a2 = ( (deltat.dt2 - deltat.dt1) / (deltat.dt1 + deltat.dt2) )
                        a3 = ( deltat.dt1 / (deltat.dt2 * (deltat.dt1 + deltat.dt2)))

                        #coef1
                        a = a2 * materialtime1.ρ[i,j-1] * mesh.vol[i,j-1]
                        b  = a2 * materialtime1.ρ[i,j] * mesh.vol[i,j]

                        num = a * dx2 + b * dx1
                        den = dx1 + dx2
                        coef1 =  num / den

                        num = velocity.v.time1[i,j-1] * dx2 + velocity.v.time1[i,j] * dx1
                        den = dx1 + dx2
                        meantime1 = num / den

                        fValue_old1 = velocity.fValues.vFaceTime1[i,j]

                        term1 = (coef1 * Df / Vf) * (-1.0) * (fValue_old1 - meantime1)

                        #coef2
                        a = a3 * materialtime2.ρ[i,j-1] * mesh.vol[i,j-1]
                        b  = a3 * materialtime2.ρ[i,j] * mesh.vol[i,j]

                        num = a * dx2 + b * dx1
                        den = dx1 + dx2
                        coef2 =  num / den

                        num = velocity.v.time2[i,j-1] * dx2 + velocity.v.time2[i,j] * dx1
                        den = dx1 + dx2
                        meantime2 = num / den

                        fValue_old2 = velocity.fValues.vFaceTime2[i,j]

                        term2 = (coef2 * Df / Vf) * (-1.0) * (fValue_old2 - meantime2)

                        #Final
                        v_rhieChow_Time[i,j] = term1 + term2

                    elseif (scheme == 3) #BDF2
                        a1 = ( 1.0 / deltat.dt1 ) + ( 1.0 / (deltat.dt1 + deltat.dt2))
                        a2 = (( 1.0 / deltat.dt1 ) + ( 1.0 / deltat.dt2))
                        a3 = -1.0 * ( deltat.dt1 / (deltat.dt2 * (deltat.dt1 + deltat.dt2)))

                        #coef1
                        a = a2 * materialtime1.ρ[i,j-1] * mesh.vol[i,j-1]
                        b  = a2 * materialtime1.ρ[i,j] * mesh.vol[i,j]

                        num = a * dx2 + b * dx1
                        den = dx1 + dx2
                        coef1 =  num / den

                        num = velocity.v.time1[i,j-1] * dx2 + velocity.v.time1[i,j] * dx1
                        den = dx1 + dx2
                        meantime1 = num / den

                        fValue_old1 = velocity.fValues.vFaceTime1[i,j]

                        term1 = (coef1 * Df / Vf) * (-1.0) * (fValue_old1 - meantime1)

                        #coef2
                        a = a3 * materialtime2.ρ[i,j-1] * mesh.vol[i,j-1]
                        b  = a3 * materialtime2.ρ[i,j] * mesh.vol[i,j]

                        num = a * dx2 + b * dx1
                        den = dx1 + dx2
                        coef2 =  num / den

                        num = velocity.v.time2[i,j-1] * dx2 + velocity.v.time2[i,j] * dx1
                        den = dx1 + dx2
                        meantime2 = num / den

                        fValue_old2 = velocity.fValues.vFaceTime2[i,j]

                        term2 = (coef2 * Df / Vf) * (-1.0) * (fValue_old2 - meantime2)

                        #Final
                        v_rhieChow_Time[i,j] = term1 + term2

                    else
                        error("Time scheme: $(scheme) unimplemented")

                    end

                end
            end
        end

    end

    return u_rhieChow_Time, v_rhieChow_Time
end

function compute_RhieChow_Time(
    AU::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,2},
    },
    AV::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,2},
    },
    AW::Union{
    SparseVector{<:AbstractFloat,<:Signed},
    SparseMatrixCSC{<:AbstractFloat,<:Signed},
    Array{<:AbstractFloat,2},
    },
    velocity::CSVelocity3D,
    mesh::UnionCSMesh3D,
    deltat::DeltaTime,
    material::CSMaterial3D;
    materialtime1::CSMaterial3D = material,
    materialtime2::CSMaterial3D = material,
    materialtime3::CSMaterial3D = material,
    T::Type{<:AbstractFloat} = Float64,
    scheme::Signed = 1,
)
    u_rhieChow_Time = zeros(T, mesh.l1, mesh.m1)
    v_rhieChow_Time = zeros(T, mesh.l1, mesh.m1)
    w_rhieChow_Time = zeros(T, mesh.l1, mesh.m1)

    #Compute u
    if mesh.l1 != 1
        for i in 2:mesh.l1
            for j in 1:mesh.m1
                for k in 1:mesh.n1
                    if velocity.u.onoff[i-1,j,k] && velocity.u.onoff[i,j,k]
                        dx1 = 0.0
                        dx2 =  0.0
                        id1 = 0
                        id2 = 0
                        D1  = 0.0
                        D2  = 0.0
                        Df = 0.0
                        Vf = 0.0
                        coef1 = 0.0
                        coef2 = 0.0
                        coef3 = 0.0
                        a = 0.0
                        b = 0.0
                        meantime1 = 0.0
                        meantime2 = 0.0
                        meantime3 = 0.0
                        fValue_old1 = 0.0
                        fValue_old2 = 0.0
                        fValue_old3 = 0.0
                        term1 = 0.0
                        term2 = 0.0
                        term3 = 0.0

                        #parameters
                        dx1 = 0.5 * mesh.dx[i-1]
                        dx2 =  0.5 * mesh.dx[i]
                        id1 = velocity.u.gIndex[i-1,j,k]
                        id2 = velocity.u.gIndex[i,j,k]

                        D1  = mesh.vol[i-1,j,k] / AU[id1,id1]
                        D2  = mesh.vol[i,j,k] / AU[id2,id2]
                        num = D1 * dx2 + D2 * dx1
                        den = dx1 + dx2
                        Df = num / den

                        num = (0.5 * mesh.dx[i-1]) * dx2 + (0.5 * mesh.dx[i]) * dx1
                        den = dx1 + dx2
                        Vf = num / den

                        if (scheme == 1) #euler
                            num = materialtime1.ρ[i-1,j,k] * mesh.vol[i-1,j,k]
                            den = deltat.dt1
                            a = num / den

                            num = materialtime1.ρ[i,j,k] * mesh.vol[i,j,k]
                            den = deltat.dt1
                            b  = num / den

                            num = a * dx2 + b * dx1
                            den = dx1 + dx2
                            coef1 =  num / den

                            num = velocity.u.time1[i-1,j,k] * dx2 + velocity.u.time1[i,j,k] * dx1
                            den = dx1 + dx2
                            meantime1 = num / den

                            fValue_old1 = velocity.fValues.uFaceTime1[i,j,k]

                            u_rhieChow_Time[i,j,k] = (coef1 * Df / Vf) * (-1.0) * (fValue_old1 - meantime1)

                        elseif (scheme == 2) #CN
                            a1 = ( deltat.dt2 / (deltat.dt1 * (deltat.dt1 + deltat.dt2)))
                            a2 = ( (deltat.dt2 - deltat.dt1) / (deltat.dt1 + deltat.dt2) )
                            a3 = ( deltat.dt1 / (deltat.dt2 * (deltat.dt1 + deltat.dt2)))

                            #coef1
                            a = a2 * materialtime1.ρ[i-1,j,k] * mesh.vol[i-1,j,k]
                            b  = a2 * materialtime1.ρ[i,j,k] * mesh.vol[i,j,k]

                            num = a * dx2 + b * dx1
                            den = dx1 + dx2
                            coef1 =  num / den

                            num = velocity.u.time1[i-1,j,k] * dx2 + velocity.u.time1[i,j,k] * dx1
                            den = dx1 + dx2
                            meantime1 = num / den

                            fValue_old1 = velocity.fValues.uFaceTime1[i,j,k]

                            term1 = (coef1 * Df / Vf) * (-1.0) * (fValue_old1 - meantime1)

                            #coef2
                            a = a3 * materialtime2.ρ[i-1,j,k] * mesh.vol[i-1,j,k]
                            b  = a3 * materialtime2.ρ[i,j,k] * mesh.vol[i,j,k]

                            num = a * dx2 + b * dx1
                            den = dx1 + dx2
                            coef2 =  num / den

                            num = velocity.u.time2[i-1,j,k] * dx2 + velocity.u.time2[i,j,k] * dx1
                            den = dx1 + dx2
                            meantime2 = num / den

                            fValue_old2 = velocity.fValues.uFaceTime2[i,j,k]

                            term2 = (coef2 * Df / Vf) * (-1.0) * (fValue_old2 - meantime2)

                            #Final
                            u_rhieChow_Time[i,j,k] = term1 + term2

                        elseif (scheme == 3) #BDF2
                            a1 = ( 1.0 / deltat.dt1 ) + ( 1.0 / (deltat.dt1 + deltat.dt2))
                            a2 = (( 1.0 / deltat.dt1 ) + ( 1.0 / deltat.dt2))
                            a3 = -1.0 * ( deltat.dt1 / (deltat.dt2 * (deltat.dt1 + deltat.dt2)))

                            #coef1
                            a = a2 * materialtime1.ρ[i-1,j,k] * mesh.vol[i-1,j,k]
                            b  = a2 * materialtime1.ρ[i,j,k] * mesh.vol[i,j,k]

                            num = a * dx2 + b * dx1
                            den = dx1 + dx2
                            coef1 =  num / den

                            num = velocity.u.time1[i-1,j,k] * dx2 + velocity.u.time1[i,j,k] * dx1
                            den = dx1 + dx2
                            meantime1 = num / den

                            fValue_old1 = velocity.fValues.uFaceTime1[i,j,k]

                            term1 = (coef1 * Df / Vf) * (-1.0) * (fValue_old1 - meantime1)

                            #coef2
                            a = a3 * materialtime2.ρ[i-1,j,k] * mesh.vol[i-1,j,k]
                            b  = a3 * materialtime2.ρ[i,j,k] * mesh.vol[i,j,k]

                            num = a * dx2 + b * dx1
                            den = dx1 + dx2
                            coef2 =  num / den

                            num = velocity.u.time2[i-1,j,k] * dx2 + velocity.u.time2[i,j,k] * dx1
                            den = dx1 + dx2
                            meantime2 = num / den

                            fValue_old2 = velocity.fValues.uFaceTime2[i,j,k]

                            term2 = (coef2 * Df / Vf) * (-1.0) * (fValue_old2 - meantime2)

                            #Final
                            u_rhieChow_Time[i,j,k] = term1 + term2

                        else
                            error("Time scheme: $(scheme) unimplemented")

                        end

                    end
                end
            end
        end

    end

    #Compute v
    if mesh.m1 != 1
        for i in 1:mesh.l1
            for j in 2:mesh.m1
                for k in 1:mesh.n1
                    if velocity.v.onoff[i,j-1,k] && velocity.v.onoff[i,j,k]
                        dx1 = 0.0
                        dx2 =  0.0
                        id1 = 0
                        id2 = 0
                        D1  = 0.0
                        D2  = 0.0
                        Df = 0.0
                        Vf = 0.0
                        coef1 = 0.0
                        coef2 = 0.0
                        coef3 = 0.0
                        a = 0.0
                        b = 0.0
                        meantime1 = 0.0
                        meantime2 = 0.0
                        meantime3 = 0.0
                        fValue_old1 = 0.0
                        fValue_old2 = 0.0
                        fValue_old3 = 0.0
                        term1 = 0.0
                        term2 = 0.0
                        term3 = 0.0

                        #parameters
                        dx1 = 0.5 * mesh.dy[j-1]
                        dx2 =  0.5 * mesh.dy[j]
                        id1 = velocity.v.gIndex[i,j-1,k]
                        id2 = velocity.v.gIndex[i,j,k]

                        D1  = mesh.vol[i,j-1,k] / AV[id1,id1]
                        D2  = mesh.vol[i,j,k] / AV[id2,id2]
                        num = D1 * dx2 + D2 * dx1
                        den = dx1 + dx2
                        Df = num / den

                        num = (0.5 * mesh.dy[j-1]) * dx2 + (0.5 * mesh.dy[j]) * dx1
                        den = dx1 + dx2
                        Vf = num / den

                        if (scheme == 1) #euler
                            num = materialtime1.ρ[i,j-1,k] * mesh.vol[i,j-1,k]
                            den = deltat.dt1
                            a = num / den

                            num = materialtime1.ρ[i,j,k] * mesh.vol[i,j,k]
                            den = deltat.dt1
                            b  = num / den

                            num = a * dx2 + b * dx1
                            den = dx1 + dx2
                            coef1 =  num / den

                            num = velocity.v.time1[i,j-1,k] * dx2 + velocity.v.time1[i,j,k] * dx1
                            den = dx1 + dx2
                            meantime1 = num / den

                            fValue_old1 = velocity.fValues.vFaceTime1[i,j,k]

                            v_rhieChow_Time[i,j,k] = (coef1 * Df / Vf) * (-1.0) * (fValue_old1 - meantime1)

                        elseif (scheme == 2) #CN
                            a1 = ( deltat.dt2 / (deltat.dt1 * (deltat.dt1 + deltat.dt2)))
                            a2 = ( (deltat.dt2 - deltat.dt1) / (deltat.dt1 + deltat.dt2) )
                            a3 = ( deltat.dt1 / (deltat.dt2 * (deltat.dt1 + deltat.dt2)))

                            #coef1
                            a = a2 * materialtime1.ρ[i,j-1,k] * mesh.vol[i,j-1,k]
                            b  = a2 * materialtime1.ρ[i,j,k] * mesh.vol[i,j,k]

                            num = a * dx2 + b * dx1
                            den = dx1 + dx2
                            coef1 =  num / den

                            num = velocity.v.time1[i,j-1,k] * dx2 + velocity.v.time1[i,j,k] * dx1
                            den = dx1 + dx2
                            meantime1 = num / den

                            fValue_old1 = velocity.fValues.vFaceTime1[i,j,k]

                            term1 = (coef1 * Df / Vf) * (-1.0) * (fValue_old1 - meantime1)

                            #coef2
                            a = a3 * materialtime2.ρ[i,j-1,k] * mesh.vol[i,j-1,k]
                            b  = a3 * materialtime2.ρ[i,j,k] * mesh.vol[i,j,k]

                            num = a * dx2 + b * dx1
                            den = dx1 + dx2
                            coef2 =  num / den

                            num = velocity.v.time2[i,j-1,k] * dx2 + velocity.v.time2[i,j,k] * dx1
                            den = dx1 + dx2
                            meantime2 = num / den

                            fValue_old2 = velocity.fValues.vFaceTime2[i,j,k]

                            term2 = (coef2 * Df / Vf) * (-1.0) * (fValue_old2 - meantime2)

                            #Final
                            v_rhieChow_Time[i,j,k] = term1 + term2

                        elseif (scheme == 3) #BDF2
                            a1 = ( 1.0 / deltat.dt1 ) + ( 1.0 / (deltat.dt1 + deltat.dt2))
                            a2 = (( 1.0 / deltat.dt1 ) + ( 1.0 / deltat.dt2))
                            a3 = -1.0 * ( deltat.dt1 / (deltat.dt2 * (deltat.dt1 + deltat.dt2)))

                            #coef1
                            a = a2 * materialtime1.ρ[i,j-1,k] * mesh.vol[i,j-1,k]
                            b  = a2 * materialtime1.ρ[i,j,k] * mesh.vol[i,j,k]

                            num = a * dx2 + b * dx1
                            den = dx1 + dx2
                            coef1 =  num / den

                            num = velocity.v.time1[i,j-1,k] * dx2 + velocity.v.time1[i,j,k] * dx1
                            den = dx1 + dx2
                            meantime1 = num / den

                            fValue_old1 = velocity.fValues.vFaceTime1[i,j,k]

                            term1 = (coef1 * Df / Vf) * (-1.0) * (fValue_old1 - meantime1)

                            #coef2
                            a = a3 * materialtime2.ρ[i,j-1,k] * mesh.vol[i,j-1,k]
                            b  = a3 * materialtime2.ρ[i,j,k] * mesh.vol[i,j,k]

                            num = a * dx2 + b * dx1
                            den = dx1 + dx2
                            coef2 =  num / den

                            num = velocity.v.time2[i,j-1,k] * dx2 + velocity.v.time2[i,j,k] * dx1
                            den = dx1 + dx2
                            meantime2 = num / den

                            fValue_old2 = velocity.fValues.vFaceTime2[i,j,k]

                            term2 = (coef2 * Df / Vf) * (-1.0) * (fValue_old2 - meantime2)

                            #Final
                            v_rhieChow_Time[i,j,k] = term1 + term2

                        else
                            error("Time scheme: $(scheme) unimplemented")

                        end

                    end
                end
            end
        end

    end


    #Compute v
    if mesh.n1 != 1
        for i in 1:mesh.l1
            for j in 1:mesh.m1
                for k in 2:mesh.n1
                    if velocity.w.onoff[i,j,k-1] && velocity.w.onoff[i,j,k]
                        dx1 = 0.0
                        dx2 =  0.0
                        id1 = 0
                        id2 = 0
                        D1  = 0.0
                        D2  = 0.0
                        Df = 0.0
                        Vf = 0.0
                        coef1 = 0.0
                        coef2 = 0.0
                        coef3 = 0.0
                        a = 0.0
                        b = 0.0
                        meantime1 = 0.0
                        meantime2 = 0.0
                        meantime3 = 0.0
                        fValue_old1 = 0.0
                        fValue_old2 = 0.0
                        fValue_old3 = 0.0
                        term1 = 0.0
                        term2 = 0.0
                        term3 = 0.0

                        #parameters
                        dx1 = 0.5 * mesh.dz[k-1]
                        dx2 =  0.5 * mesh.dz[k]
                        id1 = velocity.w.gIndex[i,j,k-1]
                        id2 = velocity.w.gIndex[i,j,k]

                        D1  = mesh.vol[i,j,k-1] / AW[id1,id1]
                        D2  = mesh.vol[i,j,k] / AW[id2,id2]
                        num = D1 * dx2 + D2 * dx1
                        den = dx1 + dx2
                        Df = num / den

                        num = (0.5 * mesh.dz[k-1]) * dx2 + (0.5 * mesh.dz[k]) * dx1
                        den = dx1 + dx2
                        Vf = num / den

                        if (scheme == 1) #euler
                            num = materialtime1.ρ[i,j,k-1] * mesh.vol[i,j,k-1]
                            den = deltat.dt1
                            a = num / den

                            num = materialtime1.ρ[i,j,k] * mesh.vol[i,j,k]
                            den = deltat.dt1
                            b  = num / den

                            num = a * dx2 + b * dx1
                            den = dx1 + dx2
                            coef1 =  num / den

                            num = velocity.w.time1[i,j,k-1] * dx2 + velocity.w.time1[i,j,k] * dx1
                            den = dx1 + dx2
                            meantime1 = num / den

                            fValue_old1 = velocity.fValues.wFaceTime1[i,j,k]

                            w_rhieChow_Time[i,j,k] = (coef1 * Df / Vf) * (-1.0) * (fValue_old1 - meantime1)

                        elseif (scheme == 2) #CN
                            a1 = ( deltat.dt2 / (deltat.dt1 * (deltat.dt1 + deltat.dt2)))
                            a2 = ( (deltat.dt2 - deltat.dt1) / (deltat.dt1 + deltat.dt2) )
                            a3 = ( deltat.dt1 / (deltat.dt2 * (deltat.dt1 + deltat.dt2)))

                            #coef1
                            a = a2 * materialtime1.ρ[i,j,k-1] * mesh.vol[i,j,k-1]
                            b  = a2 * materialtime1.ρ[i,j,k] * mesh.vol[i,j,k]

                            num = a * dx2 + b * dx1
                            den = dx1 + dx2
                            coef1 =  num / den

                            num = velocity.w.time1[i,j,k-1] * dx2 + velocity.w.time1[i,j,k] * dx1
                            den = dx1 + dx2
                            meantime1 = num / den

                            fValue_old1 = velocity.fValues.wFaceTime1[i,j,k]

                            term1 = (coef1 * Df / Vf) * (-1.0) * (fValue_old1 - meantime1)

                            #coef2
                            a = a3 * materialtime2.ρ[i,j,k-1] * mesh.vol[i,j,k-1]
                            b  = a3 * materialtime2.ρ[i,j,k] * mesh.vol[i,j,k]

                            num = a * dx2 + b * dx1
                            den = dx1 + dx2
                            coef2 =  num / den

                            num = velocity.w.time2[i,j,k-1] * dx2 + velocity.w.time2[i,j,k] * dx1
                            den = dx1 + dx2
                            meantime2 = num / den

                            fValue_old2 = velocity.fValues.wFaceTime2[i,j,k]

                            term2 = (coef2 * Df / Vf) * (-1.0) * (fValue_old2 - meantime2)

                            #Final
                            w_rhieChow_Time[i,j,k] = term1 + term2

                        elseif (scheme == 3) #BDF2
                            a1 = ( 1.0 / deltat.dt1 ) + ( 1.0 / (deltat.dt1 + deltat.dt2))
                            a2 = (( 1.0 / deltat.dt1 ) + ( 1.0 / deltat.dt2))
                            a3 = -1.0 * ( deltat.dt1 / (deltat.dt2 * (deltat.dt1 + deltat.dt2)))

                            #coef1
                            a = a2 * materialtime1.ρ[i,j,k-1] * mesh.vol[i,j,k-1]
                            b  = a2 * materialtime1.ρ[i,j,k] * mesh.vol[i,j,k]

                            num = a * dx2 + b * dx1
                            den = dx1 + dx2
                            coef1 =  num / den

                            num = velocity.w.time1[i,j,k-1] * dx2 + velocity.w.time1[i,j,k] * dx1
                            den = dx1 + dx2
                            meantime1 = num / den

                            fValue_old1 = velocity.fValues.wFaceTime1[i,j,k]

                            term1 = (coef1 * Df / Vf) * (-1.0) * (fValue_old1 - meantime1)

                            #coef2
                            a = a3 * materialtime2.ρ[i,j,k-1] * mesh.vol[i,j,k-1]
                            b  = a3 * materialtime2.ρ[i,j,k] * mesh.vol[i,j,k]

                            num = a * dx2 + b * dx1
                            den = dx1 + dx2
                            coef2 =  num / den

                            num = velocity.w.time2[i,j,k-1] * dx2 + velocity.w.time2[i,j,k] * dx1
                            den = dx1 + dx2
                            meantime2 = num / den

                            fValue_old2 = velocity.fValues.wFaceTime2[i,j,k]

                            term2 = (coef2 * Df / Vf) * (-1.0) * (fValue_old2 - meantime2)

                            #Final
                            w_rhieChow_Time[i,j,k] = term1 + term2

                        else
                            error("Time scheme: $(scheme) unimplemented")

                        end

                    end
                end
            end
        end

    end

    return u_rhieChow_Time, v_rhieChow_Time, w_rhieChow_Time
end
