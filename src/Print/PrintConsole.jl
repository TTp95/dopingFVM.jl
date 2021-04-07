"""

"""
function print_consoleError(
    system::SystemControl,
    args...;
    mass::AbstractFloat = -1.0,
)
    nt = length(args)

    if (system.iterationsTimeStep == 1)
        print("\n")
        @printf("%8s  %6s ", "Time", "n")

        for x in 1:nt
            if (typeof(args[x]) <: UnionCSVelocity)
                if (typeof(args[x]) <: CSVelocity1D)
                    @printf(" %-16s ", "error $(args[x].u.key)")

                elseif (typeof(args[x]) <: CSVelocity2D)
                    @printf(" %-16s ", "error $(args[x].u.key)")
                    @printf(" %-16s ", "error $(args[x].v.key)")

                elseif (typeof(args[x]) <: CSVelocity3D)
                    @printf(" %-16s ", "error $(args[x].u.key)")
                    @printf(" %-16s ", "error $(args[x].v.key)")
                    @printf(" %-16s ", "error $(args[x].w.key)")

                end

            elseif (typeof(args[x]) <: PhiCartesianStructured)
                @printf(" %-16s ", "error $(args[x].key)")

            else
                error("Argumment number $(x + 1) incompatible...")
            end
        end

        if (mass >= 0.0)
            @printf(" %-16s ", "SMAX")
        end

        @printf("\n")
        @printf("\n")
    end

    @printf(" %8.4g  %6d ", system.time, system.iterationsTimeStep)

    for x in 1:nt
        if (typeof(args[x]) <: UnionCSVelocity)
            if (typeof(args[x]) <: CSVelocity1D)
                @printf(" %-16.8g ", args[x].u.errorI)

            elseif (typeof(args[x]) <: CSVelocity2D)
                @printf(" %-16.8g ", args[x].u.errorI)
                @printf(" %-16.8g ", args[x].v.errorI)

            elseif (typeof(args[x]) <: CSVelocity3D)
                @printf(" %-16.8g ", args[x].u.errorI)
                @printf(" %-16.8g ", args[x].v.errorI)
                @printf(" %-16.8g ", args[x].w.errorI)

            end

        elseif (typeof(args[x]) <: PhiCartesianStructured)
            @printf(" %-16.8g ", args[x].errorI)

        else
            error("Argumment number $(x + 1) incompatible...")
        end
    end

    if (mass >= 0.0)
        @printf(" %-16.8g ", mass)
    end

    @printf("\n")


    return nothing
end
