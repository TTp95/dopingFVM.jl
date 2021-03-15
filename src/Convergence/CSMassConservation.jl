"""

"""
function mass_conservation end

function mass_conservationCS(
    velocity::CSVelocity1D,
    material::CSMaterial1D,
    mesh::UnionCSMesh1D;
    threads = false,
)
    #Inicializate variables
    massFlux = zeros(mesh.l1)

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            xflux = 0.0

            #Flux in x-cord
            if (i == 1)
                dens1 = material.ρ[i]
                dens2 = density_interpolationCS([mesh.dx[i]; mesh.dx[i+1]], [material.ρ[i]; material.ρ[i+1]])
                xflux1 = velocity.fValues.uFace[i] * dens1
                xflux2 = velocity.fValues.uFace[i+1] * dens2
                xflux = xflux2 - xflux1
            elseif (i == mesh.l1)
                dens1 = density_interpolationCS([mesh.dx[i]; mesh.dx[i-1]], [material.ρ[i]; material.ρ[i-1]])
                dens2 = material.ρ[i]
                xflux1 = velocity.fValues.uFace[i] * dens1
                xflux2 = velocity.fValues.uFace[i+1] * dens2
                xflux = xflux2 - xflux1
            else
                dens1 = density_interpolationCS([mesh.dx[i]; mesh.dx[i-1]], [material.ρ[i]; material.ρ[i-1]])
                dens2 = density_interpolationCS([mesh.dx[i]; mesh.dx[i+1]], [material.ρ[i]; material.ρ[i+1]])
                xflux1 = velocity.fValues.uFace[i] * dens1
                xflux2 = velocity.fValues.uFace[i+1] * dens2
                xflux = xflux2 - xflux1
            end

            @inbounds massFlux[i] = abs(xflux + yflux)
        end
    elseif !threads
        for i in 1:mesh.l1
            xflux = 0.0

            #Flux in x-cord
            if (i == 1)
                dens1 = material.ρ[i]
                dens2 = density_interpolationCS([mesh.dx[i]; mesh.dx[i+1]], [material.ρ[i]; material.ρ[i+1]])
                xflux1 = velocity.fValues.uFace[i] * dens1
                xflux2 = velocity.fValues.uFace[i+1] * dens2
                xflux = xflux2 - xflux1
            elseif (i == mesh.l1)
                dens1 = density_interpolationCS([mesh.dx[i]; mesh.dx[i-1]], [material.ρ[i]; material.ρ[i-1]])
                dens2 = material.ρ[i]
                xflux1 = velocity.fValues.uFace[i] * dens1
                xflux2 = velocity.fValues.uFace[i+1] * dens2
                xflux = xflux2 - xflux1
            else
                dens1 = density_interpolationCS([mesh.dx[i]; mesh.dx[i-1]], [material.ρ[i]; material.ρ[i-1]])
                dens2 = density_interpolationCS([mesh.dx[i]; mesh.dx[i+1]], [material.ρ[i]; material.ρ[i+1]])
                xflux1 = velocity.fValues.uFace[i] * dens1
                xflux2 = velocity.fValues.uFace[i+1] * dens2
                xflux = xflux2 - xflux1
            end

            @inbounds massFlux[i] = abs(xflux + yflux)
        end
    end

    return maximum(massFlux)
end

function mass_conservationCS(
    velocity::CSVelocity2D,
    material::CSMaterial2D,
    mesh::UnionCSMesh2D;
    threads = false,
)
    #Inicializate variables
    massFlux = zeros(mesh.l1, mesh.m1)

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            for j in 1:mesh.m1
                xflux = 0.0
                yflux = 0.0

                #Flux in x-cord
                if (i == 1)
                    dens1 = material.ρ[i,j]
                    dens2 = density_interpolationCS([mesh.dx[i]; mesh.dx[i+1]], [material.ρ[i,j]; material.ρ[i+1,j]])
                    xflux1 = velocity.fValues.uFace[i,j] * dens1 * mesh.dy[j]
                    xflux2 = velocity.fValues.uFace[i+1,j] * dens2 * mesh.dy[j]
                    xflux = xflux2 - xflux1
                elseif (i == mesh.l1)
                    dens1 = density_interpolationCS([mesh.dx[i]; mesh.dx[i-1]], [material.ρ[i,j]; material.ρ[i-1,j]])
                    dens2 = material.ρ[i,j]
                    xflux1 = velocity.fValues.uFace[i,j] * dens1 * mesh.dy[j]
                    xflux2 = velocity.fValues.uFace[i+1,j] * dens2 * mesh.dy[j]
                    xflux = xflux2 - xflux1
                else
                    dens1 = density_interpolationCS([mesh.dx[i]; mesh.dx[i-1]], [material.ρ[i,j]; material.ρ[i-1,j]])
                    dens2 = density_interpolationCS([mesh.dx[i]; mesh.dx[i+1]], [material.ρ[i,j]; material.ρ[i+1,j]])
                    xflux1 = velocity.fValues.uFace[i,j] * dens1 * mesh.dy[j]
                    xflux2 = velocity.fValues.uFace[i+1,j] * dens2 * mesh.dy[j]
                    xflux = xflux2 - xflux1
                end

                #Flux in y-cord
                if (j == 1)
                    dens1 = material.ρ[i,j]
                    dens2 = density_interpolationCS([mesh.dy[j]; mesh.dy[j+1]], [material.ρ[i,j]; material.ρ[i,j+1]])
                    yflux1 = velocity.fValues.vFace[i,j] * dens1 * mesh.dx[i]
                    yflux2 = velocity.fValues.vFace[i,j+1] * dens2 * mesh.dx[i]
                    yflux = yflux2 - yflux1
                elseif (j == mesh.m1)
                    dens1 = density_interpolationCS([mesh.dy[j]; mesh.dy[j-1]], [material.ρ[i,j]; material.ρ[i,j-1]])
                    dens2 = material.ρ[i,j]
                    yflux1 = velocity.fValues.vFace[i,j,k] * dens1 * mesh.dx[i]
                    yflux2 = velocity.fValues.vFace[i,j+1] * dens2 * mesh.dx[i]
                    yflux = yflux2 - yflux1
                else
                    dens1 = density_interpolationCS([mesh.dy[j]; mesh.dy[j-1]], [material.ρ[i,j]; material.ρ[i,j-1]])
                    dens2 = density_interpolationCS([mesh.dy[j]; mesh.dy[j+1]], [material.ρ[i,j]; material.ρ[i,j+1]])
                    yflux1 = velocity.fValues.vFace[i,j] * dens1 * mesh.dx[i]
                    yflux2 = velocity.fValues.vFace[i,j+1] * dens2 * mesh.dx[i]
                    yflux = yflux2 - yflux1
                end

                @inbounds massFlux[i,j] = abs(xflux + yflux)
            end
        end
    elseif !threads
        for i in 1:mesh.l1
            for j in 1:mesh.m1
                xflux = 0.0
                yflux = 0.0

                #Flux in x-cord
                if (i == 1)
                    dens1 = material.ρ[i,j]
                    dens2 = density_interpolationCS([mesh.dx[i]; mesh.dx[i+1]], [material.ρ[i,j]; material.ρ[i+1,j]])
                    xflux1 = velocity.fValues.uFace[i,j] * dens1 * mesh.dy[j]
                    xflux2 = velocity.fValues.uFace[i+1,j] * dens2 * mesh.dy[j]
                    xflux = xflux2 - xflux1
                elseif (i == mesh.l1)
                    dens1 = density_interpolationCS([mesh.dx[i]; mesh.dx[i-1]], [material.ρ[i,j]; material.ρ[i-1,j]])
                    dens2 = material.ρ[i,j]
                    xflux1 = velocity.fValues.uFace[i,j] * dens1 * mesh.dy[j]
                    xflux2 = velocity.fValues.uFace[i+1,j] * dens2 * mesh.dy[j]
                    xflux = xflux2 - xflux1
                else
                    dens1 = density_interpolationCS([mesh.dx[i]; mesh.dx[i-1]], [material.ρ[i,j]; material.ρ[i-1,j]])
                    dens2 = density_interpolationCS([mesh.dx[i]; mesh.dx[i+1]], [material.ρ[i,j]; material.ρ[i+1,j]])
                    xflux1 = velocity.fValues.uFace[i,j] * dens1 * mesh.dy[j]
                    xflux2 = velocity.fValues.uFace[i+1,j] * dens2 * mesh.dy[j]
                    xflux = xflux2 - xflux1
                end

                #Flux in y-cord
                if (j == 1)
                    dens1 = material.ρ[i,j]
                    dens2 = density_interpolationCS([mesh.dy[j]; mesh.dy[j+1]], [material.ρ[i,j]; material.ρ[i,j+1]])
                    yflux1 = velocity.fValues.vFace[i,j] * dens1 * mesh.dx[i]
                    yflux2 = velocity.fValues.vFace[i,j+1] * dens2 * mesh.dx[i]
                    yflux = yflux2 - yflux1
                elseif (j == mesh.m1)
                    dens1 = density_interpolationCS([mesh.dy[j]; mesh.dy[j-1]], [material.ρ[i,j]; material.ρ[i,j-1]])
                    dens2 = material.ρ[i,j]
                    yflux1 = velocity.fValues.vFace[i,j,k] * dens1 * mesh.dx[i]
                    yflux2 = velocity.fValues.vFace[i,j+1] * dens2 * mesh.dx[i]
                    yflux = yflux2 - yflux1
                else
                    dens1 = density_interpolationCS([mesh.dy[j]; mesh.dy[j-1]], [material.ρ[i,j]; material.ρ[i,j-1]])
                    dens2 = density_interpolationCS([mesh.dy[j]; mesh.dy[j+1]], [material.ρ[i,j]; material.ρ[i,j+1]])
                    yflux1 = velocity.fValues.vFace[i,j] * dens1 * mesh.dx[i]
                    yflux2 = velocity.fValues.vFace[i,j+1] * dens2 * mesh.dx[i]
                    yflux = yflux2 - yflux1
                end

                @inbounds massFlux[i,j] = abs(xflux + yflux)
            end
        end
    end

    return maximum(massFlux)
end

function mass_conservationCS(
    velocity::CSVelocity3D,
    material::CSMaterial3D,
    mesh::UnionCSMesh3D;
    threads = false,
)
    #Inicializate variables
    massFlux = zeros(mesh.l1, mesh.m1, mesh.n1)

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            for j in 1:mesh.m1
                for k in 1:mesh.n1

                    xflux = 0.0
                    yflux = 0.0
                    zflux = 0.0

                    #Flux in x-cord
                    if (i == 1)
                        dens1 = material.ρ[i,j,k]
                        dens2 = density_interpolationCS([mesh.dx[i]; mesh.dx[i+1]], [material.ρ[i,j,k]; material.ρ[i+1,j,k]])
                        xflux1 = velocity.fValues.uFace[i,j,k] * dens1 * mesh.dy[j] * mesh.dz[k]
                        xflux2 = velocity.fValues.uFace[i+1,j,k] * dens2 * mesh.dy[j] * mesh.dz[k]
                        xflux = xflux2 - xflux1
                    elseif (i == mesh.l1)
                        dens1 = density_interpolationCS([mesh.dx[i]; mesh.dx[i-1]], [material.ρ[i,j,k]; material.ρ[i-1,j,k]])
                        dens2 = material.ρ[i,j,k]
                        xflux1 = velocity.fValues.uFace[i,j,k] * dens1 * mesh.dy[j] * mesh.dz[k]
                        xflux2 = velocity.fValues.uFace[i+1,j,k] * dens2 * mesh.dy[j] * mesh.dz[k]
                        xflux = xflux2 - xflux1
                    else
                        dens1 = density_interpolationCS([mesh.dx[i]; mesh.dx[i-1]], [material.ρ[i,j,k]; material.ρ[i-1,j,k]])
                        dens2 = density_interpolationCS([mesh.dx[i]; mesh.dx[i+1]], [material.ρ[i,j,k]; material.ρ[i+1,j,k]])
                        xflux1 = velocity.fValues.uFace[i,j,k] * dens1 * mesh.dy[j] * mesh.dz[k]
                        xflux2 = velocity.fValues.uFace[i+1,j,k] * dens2 * mesh.dy[j] * mesh.dz[k]
                        xflux = xflux2 - xflux1
                    end

                    #Flux in y-cord
                    if (j == 1)
                        dens1 = material.ρ[i,j,k]
                        dens2 = density_interpolationCS([mesh.dy[j]; mesh.dy[j+1]], [material.ρ[i,j,k]; material.ρ[i,j+1,k]])
                        yflux1 = velocity.fValues.vFace[i,j,k] * dens1 * mesh.dx[i] * mesh.dz[k]
                        yflux2 = velocity.fValues.vFace[i,j+1,k] * dens2 * mesh.dx[i] * mesh.dz[k]
                        yflux = yflux2 - yflux1
                    elseif (j == mesh.m1)
                        dens1 = density_interpolationCS([mesh.dy[j]; mesh.dy[j-1]], [material.ρ[i,j,k]; material.ρ[i,j-1,k]])
                        dens2 = material.ρ[i,j,k]
                        yflux1 = velocity.fValues.vFace[i,j,k] * dens1 * mesh.dx[i] * mesh.dz[k]
                        yflux2 = velocity.fValues.vFace[i,j+1,k] * dens2 * mesh.dx[i] * mesh.dz[k]
                        yflux = yflux2 - yflux1
                    else
                        dens1 = density_interpolationCS([mesh.dy[j]; mesh.dy[j-1]], [material.ρ[i,j,k]; material.ρ[i,j-1,k]])
                        dens2 = density_interpolationCS([mesh.dy[j]; mesh.dy[j+1]], [material.ρ[i,j,k]; material.ρ[i,j+1,k]])
                        yflux1 = velocity.fValues.vFace[i,j,k] * dens1 * mesh.dx[i] * mesh.dz[k]
                        yflux2 = velocity.fValues.vFace[i,j+1,k] * dens2 * mesh.dx[i] * mesh.dz[k]
                        yflux = yflux2 - yflux1
                    end

                    #Flux in z-cord
                    if (k == 1)
                        dens1 = material.ρ[i,j,k]
                        dens2 = density_interpolationCS([mesh.dz[k]; mesh.dz[k+1]], [material.ρ[i,j,k]; material.ρ[i,j,k+1]])
                        zflux1 = velocity.fValues.wFace[i,j,k] * dens1 * mesh.dx[i] * mesh.dy[j]
                        zflux2 = velocity.fValues.wFace[i,j,k+1] * dens2 * mesh.dx[i] * mesh.dy[j]
                        zflux = zflux2 - zflux1
                    elseif (k == mesh.n1)
                        dens1 = density_interpolationCS([mesh.dz[k]; mesh.dz[k-1]], [material.ρ[i,j,k]; material.ρ[i,j,k-1]])
                        dens2 = material.ρ[i,j,k]
                        zflux1 = velocity.fValues.wFace[i,j,k] * dens1 * mesh.dx[i] * mesh.dy[j]
                        zflux2 = velocity.fValues.wFace[i,j,k+1] * dens2 * mesh.dx[i] * mesh.dy[j]
                        zflux = zflux2 - zflux1
                    else
                        dens1 = density_interpolationCS([mesh.dz[k]; mesh.dz[k-1]], [material.ρ[i,j,k]; material.ρ[i,j,k-1]])
                        dens2 = density_interpolationCS([mesh.dz[k]; mesh.dz[k+1]], [material.ρ[i,j,k]; material.ρ[i,j,k+1]])
                        zflux1 = velocity.fValues.wFace[i,j,k] * dens1 * mesh.dx[i] * mesh.dy[j]
                        zflux2 = velocity.fValues.wFace[i,j,k+1] * dens2 * mesh.dx[i] * mesh.dy[j]
                        zflux = zflux2 - zflux1
                    end

                    @inbounds massFlux[i,j,k] = abs(xflux + yflux + zflux)
                end
            end
        end
    elseif !threads
        for i in 1:mesh.l1
            for j in 1:mesh.m1
                for k in 1:mesh.n1

                    xflux = 0.0
                    yflux = 0.0
                    zflux = 0.0

                    #Flux in x-cord
                    if (i == 1)
                        dens1 = material.ρ[i,j,k]
                        dens2 = density_interpolationCS([mesh.dx[i]; mesh.dx[i+1]], [material.ρ[i,j,k]; material.ρ[i+1,j,k]])
                        xflux1 = velocity.fValues.uFace[i,j,k] * dens1 * mesh.dy[j] * mesh.dz[k]
                        xflux2 = velocity.fValues.uFace[i+1,j,k] * dens2 * mesh.dy[j] * mesh.dz[k]
                        xflux = xflux2 - xflux1
                    elseif (i == mesh.l1)
                        dens1 = density_interpolationCS([mesh.dx[i]; mesh.dx[i-1]], [material.ρ[i,j,k]; material.ρ[i-1,j,k]])
                        dens2 = material.ρ[i,j,k]
                        xflux1 = velocity.fValues.uFace[i,j,k] * dens1 * mesh.dy[j] * mesh.dz[k]
                        xflux2 = velocity.fValues.uFace[i+1,j,k] * dens2 * mesh.dy[j] * mesh.dz[k]
                        xflux = xflux2 - xflux1
                    else
                        dens1 = density_interpolationCS([mesh.dx[i]; mesh.dx[i-1]], [material.ρ[i,j,k]; material.ρ[i-1,j,k]])
                        dens2 = density_interpolationCS([mesh.dx[i]; mesh.dx[i+1]], [material.ρ[i,j,k]; material.ρ[i+1,j,k]])
                        xflux1 = velocity.fValues.uFace[i,j,k] * dens1 * mesh.dy[j] * mesh.dz[k]
                        xflux2 = velocity.fValues.uFace[i+1,j,k] * dens2 * mesh.dy[j] * mesh.dz[k]
                        xflux = xflux2 - xflux1
                    end

                    #Flux in y-cord
                    if (j == 1)
                        dens1 = material.ρ[i,j,k]
                        dens2 = density_interpolationCS([mesh.dy[j]; mesh.dy[j+1]], [material.ρ[i,j,k]; material.ρ[i,j+1,k]])
                        yflux1 = velocity.fValues.vFace[i,j,k] * dens1 * mesh.dx[i] * mesh.dz[k]
                        yflux2 = velocity.fValues.vFace[i,j+1,k] * dens2 * mesh.dx[i] * mesh.dz[k]
                        yflux = yflux2 - yflux1
                    elseif (j == mesh.m1)
                        dens1 = density_interpolationCS([mesh.dy[j]; mesh.dy[j-1]], [material.ρ[i,j,k]; material.ρ[i,j-1,k]])
                        dens2 = material.ρ[i,j,k]
                        yflux1 = velocity.fValues.vFace[i,j,k] * dens1 * mesh.dx[i] * mesh.dz[k]
                        yflux2 = velocity.fValues.vFace[i,j+1,k] * dens2 * mesh.dx[i] * mesh.dz[k]
                        yflux = yflux2 - yflux1
                    else
                        dens1 = density_interpolationCS([mesh.dy[j]; mesh.dy[j-1]], [material.ρ[i,j,k]; material.ρ[i,j-1,k]])
                        dens2 = density_interpolationCS([mesh.dy[j]; mesh.dy[j+1]], [material.ρ[i,j,k]; material.ρ[i,j+1,k]])
                        yflux1 = velocity.fValues.vFace[i,j,k] * dens1 * mesh.dx[i] * mesh.dz[k]
                        yflux2 = velocity.fValues.vFace[i,j+1,k] * dens2 * mesh.dx[i] * mesh.dz[k]
                        yflux = yflux2 - yflux1
                    end

                    #Flux in z-cord
                    if (k == 1)
                        dens1 = material.ρ[i,j,k]
                        dens2 = density_interpolationCS([mesh.dz[k]; mesh.dz[k+1]], [material.ρ[i,j,k]; material.ρ[i,j,k+1]])
                        zflux1 = velocity.fValues.wFace[i,j,k] * dens1 * mesh.dx[i] * mesh.dy[j]
                        zflux2 = velocity.fValues.wFace[i,j,k+1] * dens2 * mesh.dx[i] * mesh.dy[j]
                        zflux = zflux2 - zflux1
                    elseif (k == mesh.n1)
                        dens1 = density_interpolationCS([mesh.dz[k]; mesh.dz[k-1]], [material.ρ[i,j,k]; material.ρ[i,j,k-1]])
                        dens2 = material.ρ[i,j,k]
                        zflux1 = velocity.fValues.wFace[i,j,k] * dens1 * mesh.dx[i] * mesh.dy[j]
                        zflux2 = velocity.fValues.wFace[i,j,k+1] * dens2 * mesh.dx[i] * mesh.dy[j]
                        zflux = zflux2 - zflux1
                    else
                        dens1 = density_interpolationCS([mesh.dz[k]; mesh.dz[k-1]], [material.ρ[i,j,k]; material.ρ[i,j,k-1]])
                        dens2 = density_interpolationCS([mesh.dz[k]; mesh.dz[k+1]], [material.ρ[i,j,k]; material.ρ[i,j,k+1]])
                        zflux1 = velocity.fValues.wFace[i,j,k] * dens1 * mesh.dx[i] * mesh.dy[j]
                        zflux2 = velocity.fValues.wFace[i,j,k+1] * dens2 * mesh.dx[i] * mesh.dy[j]
                        zflux = zflux2 - zflux1
                    end

                    @inbounds massFlux[i,j,k] = abs(xflux + yflux + zflux)
                end
            end
        end
    end

    return maximum(massFlux)
end

function mass_conservationCS(
    velocity::CSVelocity1D,
    material::UnionConstantMaterial,
    mesh::UnionCSMesh1D;
    threads = false,
)
    #Inicializate variables
    massFlux = zeros(mesh.l1)

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            xflux = 0.0

            #Flux in x-cord
            dens1 = material.ρ
            dens2 = material.ρ
            xflux1 = velocity.fValues.uFace[i] * dens1
            xflux2 = velocity.fValues.uFace[i+1] * dens2
            xflux = xflux2 - xflux1

            @inbounds massFlux[i] = abs(xflux + yflux)
        end
    elseif !threads
        for i in 1:mesh.l1
            xflux = 0.0

            #Flux in x-cord
            dens1 = material.ρ
            dens2 = material.ρ
            xflux1 = velocity.fValues.uFace[i] * dens1
            xflux2 = velocity.fValues.uFace[i+1] * dens2
            xflux = xflux2 - xflux1

            @inbounds massFlux[i] = abs(xflux + yflux)
        end
    end

    return maximum(massFlux)
end

function mass_conservationCS(
    velocity::CSVelocity2D,
    material::UnionConstantMaterial,
    mesh::UnionCSMesh2D;
    threads = false,
)
    #Inicializate variables
    massFlux = zeros(mesh.l1, mesh.m1)

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            for j in 1:mesh.m1
                xflux = 0.0
                yflux = 0.0

                #Flux in x-cord
                dens1 = material.ρ
                dens2 = material.ρ
                xflux1 = velocity.fValues.uFace[i,j] * dens1 * mesh.dy[j]
                xflux2 = velocity.fValues.uFace[i+1,j] * dens2 * mesh.dy[j]
                xflux = xflux2 - xflux1

                #Flux in y-cord
                dens1 = material.ρ
                dens2 = material.ρ
                yflux1 = velocity.fValues.vFace[i,j] * dens1 * mesh.dx[i]
                yflux2 = velocity.fValues.vFace[i,j+1] * dens2 * mesh.dx[i]
                yflux = yflux2 - yflux1

                @inbounds massFlux[i,j] = abs(xflux + yflux)
            end
        end
    elseif !threads
        for i in 1:mesh.l1
            for j in 1:mesh.m1
                xflux = 0.0
                yflux = 0.0

                #Flux in x-cord
                dens1 = material.ρ
                dens2 = material.ρ
                xflux1 = velocity.fValues.uFace[i,j] * dens1 * mesh.dy[j]
                xflux2 = velocity.fValues.uFace[i+1,j] * dens2 * mesh.dy[j]
                xflux = xflux2 - xflux1

                #Flux in y-cord
                dens1 = material.ρ
                dens2 = material.ρ
                yflux1 = velocity.fValues.vFace[i,j] * dens1 * mesh.dx[i]
                yflux2 = velocity.fValues.vFace[i,j+1] * dens2 * mesh.dx[i]
                yflux = yflux2 - yflux1

                @inbounds massFlux[i,j] = abs(xflux + yflux)
            end
        end
    end

    return maximum(massFlux)
end

function mass_conservationCS(
    velocity::CSVelocity3D,
    material::UnionConstantMaterial,
    mesh::UnionCSMesh3D;
    threads = false,
)
    #Inicializate variables
    massFlux = zeros(mesh.l1, mesh.m1, mesh.n1)

    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            for j in 1:mesh.m1
                for k in 1:mesh.n1

                    xflux = 0.0
                    yflux = 0.0
                    zflux = 0.0

                    #Flux in x-cord
                    dens1 = material.ρ
                    dens2 = material.ρ
                    xflux1 = velocity.fValues.uFace[i,j,k] * dens1 * mesh.dy[j] * mesh.dz[k]
                    xflux2 = velocity.fValues.uFace[i+1,j,k] * dens2 * mesh.dy[j] * mesh.dz[k]
                    xflux = xflux2 - xflux1

                    #Flux in y-cord
                    dens1 = material.ρ
                    dens2 = material.ρ
                    yflux1 = velocity.fValues.vFace[i,j,k] * dens1 * mesh.dx[i] * mesh.dz[k]
                    yflux2 = velocity.fValues.vFace[i,j+1,k] * dens2 * mesh.dx[i] * mesh.dz[k]
                    yflux = yflux2 - yflux1

                    #Flux in z-cord
                    dens1 = material.ρ
                    dens2 = material.ρ
                    zflux1 = velocity.fValues.wFace[i,j,k] * dens1 * mesh.dx[i] * mesh.dy[j]
                    zflux2 = velocity.fValues.wFace[i,j,k+1] * dens2 * mesh.dx[i] * mesh.dy[j]
                    zflux = zflux2 - zflux1

                    @inbounds massFlux[i,j,k] = abs(xflux + yflux + zflux)
                end
            end
        end
    elseif !threads
        for i in 1:mesh.l1
            for j in 1:mesh.m1
                for k in 1:mesh.n1

                    xflux = 0.0
                    yflux = 0.0
                    zflux = 0.0

                    #Flux in x-cord
                    dens1 = material.ρ
                    dens2 = material.ρ
                    xflux1 = velocity.fValues.uFace[i,j,k] * dens1 * mesh.dy[j] * mesh.dz[k]
                    xflux2 = velocity.fValues.uFace[i+1,j,k] * dens2 * mesh.dy[j] * mesh.dz[k]
                    xflux = xflux2 - xflux1

                    #Flux in y-cord
                    dens1 = material.ρ
                    dens2 = material.ρ
                    yflux1 = velocity.fValues.vFace[i,j,k] * dens1 * mesh.dx[i] * mesh.dz[k]
                    yflux2 = velocity.fValues.vFace[i,j+1,k] * dens2 * mesh.dx[i] * mesh.dz[k]
                    yflux = yflux2 - yflux1

                    #Flux in z-cord
                    dens1 = material.ρ
                    dens2 = material.ρ
                    zflux1 = velocity.fValues.wFace[i,j,k] * dens1 * mesh.dx[i] * mesh.dy[j]
                    zflux2 = velocity.fValues.wFace[i,j,k+1] * dens2 * mesh.dx[i] * mesh.dy[j]
                    zflux = zflux2 - zflux1

                    @inbounds massFlux[i,j,k] = abs(xflux + yflux + zflux)
                end
            end
        end
    end

    return maximum(massFlux)
end
