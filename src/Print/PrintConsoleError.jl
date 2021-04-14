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
        print("Time  n")

        for x in 1:nt
            if (typeof(args[x]) <: UnionCSVelocity)
                if (typeof(args[x]) <: CSVelocity1D)
                    print(" error $(args[x].u.key) ")

                elseif (typeof(args[x]) <: CSVelocity2D)
                    print(" error $(args[x].u.key) ")
                    print(" error $(args[x].v.key) ")

                elseif (typeof(args[x]) <: CSVelocity3D)
                    print(" error $(args[x].u.key) ")
                    print(" error $(args[x].v.key) ")
                    print(" error $(args[x].w.key) ")

                end

            elseif (typeof(args[x]) <: PhiCartesianStructured)
                print(" error $(args[x].key) ")

            else
                error("Argumment number $(x + 1) incompatible...")
            end
        end

        if (mass >= 0.0)
            print(" SMAX ")
        end

        print("\n")
        print("\n")
    end

    print("$(system.time)  $(system.iterationsTimeStep)")

    for x in 1:nt
        if (typeof(args[x]) <: UnionCSVelocity)
            if (typeof(args[x]) <: CSVelocity1D)
                print(" $(args[x].u.errorI) ")

            elseif (typeof(args[x]) <: CSVelocity2D)
                print(" $(args[x].u.errorI) ")
                print(" $(args[x].v.errorI) ")

            elseif (typeof(args[x]) <: CSVelocity3D)
                print(" $(args[x].u.errorI) ")
                print(" $(args[x].v.errorI) ")
                print(" $(args[x].w.errorI) ")

            end

        elseif (typeof(args[x]) <: PhiCartesianStructured)
            print(" $(args[x].errorI) ")

        else
            error("Argumment number $(x + 1) incompatible...")
        end
    end

    if (mass >= 0.0)
        print(" $mass ")
    end

    print("\n")


    return nothing
end
