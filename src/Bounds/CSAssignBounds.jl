"""

"""
function assign_bounds! end

function assign_bounds!(
    phi::CSPhi1D,
    mesh::UnionCSMesh1D,
    west::BoundsStructured,
    east::BoundsStructured,
    dict::Dict{String,BoundsStructured};
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
    threads = false,
)
    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            if phi.bounds[i]
                for nn in 1:phi.nbounds[i]
                    gid = phi.gIndex[i]
                    if (dict["$(gid)g$(nn)"].cord == 'w')
                        dict["$(gid)g$(nn)"].kind = west.kind
                        dict["$(gid)g$(nn)"].cord = west.cord
                        dict["$(gid)g$(nn)"].ρ = west.ρ
                        dict["$(gid)g$(nn)"].Γ = west.Γ
                        dict["$(gid)g$(nn)"].a = west.a
                        dict["$(gid)g$(nn)"].b = west.b
                        dict["$(gid)g$(nn)"].c = west.c
                        dict["$(gid)g$(nn)"].eval = west.eval
                    end

                    if (dict["$(gid)g$(nn)"].cord == 'e')
                        dict["$(gid)g$(nn)"].kind = east.kind
                        dict["$(gid)g$(nn)"].cord = east.cord
                        dict["$(gid)g$(nn)"].ρ = east.ρ
                        dict["$(gid)g$(nn)"].Γ = east.Γ
                        dict["$(gid)g$(nn)"].a = east.a
                        dict["$(gid)g$(nn)"].b = east.b
                        dict["$(gid)g$(nn)"].c = east.c
                        dict["$(gid)g$(nn)"].eval = east.eval
                    end
                end
            end
        end
    elseif !threads
        for i in 1:mesh.l1
            if phi.bounds[i]
                for nn in 1:phi.nbounds[i]
                    gid = phi.gIndex[i]
                    if (dict["$(gid)g$(nn)"].cord == 'w')
                        dict["$(gid)g$(nn)"].kind = west.kind
                        dict["$(gid)g$(nn)"].cord = west.cord
                        dict["$(gid)g$(nn)"].ρ = west.ρ
                        dict["$(gid)g$(nn)"].Γ = west.Γ
                        dict["$(gid)g$(nn)"].a = west.a
                        dict["$(gid)g$(nn)"].b = west.b
                        dict["$(gid)g$(nn)"].c = west.c
                        dict["$(gid)g$(nn)"].eval = west.eval
                    end

                    if (dict["$(gid)g$(nn)"].cord == 'e')
                        dict["$(gid)g$(nn)"].kind = east.kind
                        dict["$(gid)g$(nn)"].cord = east.cord
                        dict["$(gid)g$(nn)"].ρ = east.ρ
                        dict["$(gid)g$(nn)"].Γ = east.Γ
                        dict["$(gid)g$(nn)"].a = east.a
                        dict["$(gid)g$(nn)"].b = east.b
                        dict["$(gid)g$(nn)"].c = east.c
                        dict["$(gid)g$(nn)"].eval = east.eval
                    end
                end
            end
        end
    end

    return dict
end

function assign_bounds!(
    phi::CSPhi2D,
    mesh::UnionCSMesh2D,
    west::BoundsStructured,
    east::BoundsStructured,
    south::BoundsStructured,
    north::BoundsStructured,
    dict::Dict{String,BoundsStructured};
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
    threads = false,
)
    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            for j in 1:mesh.m1
                if phi.bounds[i,j]
                    for nn in 1:phi.nbounds[i,j]
                        gid = phi.gIndex[i,j]
                        if (dict["$(gid)g$(nn)"].cord == 'w')
                            dict["$(gid)g$(nn)"].kind = west.kind
                            dict["$(gid)g$(nn)"].cord = west.cord
                            dict["$(gid)g$(nn)"].ρ = west.ρ
                            dict["$(gid)g$(nn)"].Γ = west.Γ
                            dict["$(gid)g$(nn)"].a = west.a
                            dict["$(gid)g$(nn)"].b = west.b
                            dict["$(gid)g$(nn)"].c = west.c
                            dict["$(gid)g$(nn)"].eval = west.eval
                        end

                        if (dict["$(gid)g$(nn)"].cord == 'e')
                            dict["$(gid)g$(nn)"].kind = east.kind
                            dict["$(gid)g$(nn)"].cord = east.cord
                            dict["$(gid)g$(nn)"].ρ = east.ρ
                            dict["$(gid)g$(nn)"].Γ = east.Γ
                            dict["$(gid)g$(nn)"].a = east.a
                            dict["$(gid)g$(nn)"].b = east.b
                            dict["$(gid)g$(nn)"].c = east.c
                            dict["$(gid)g$(nn)"].eval = east.eval
                        end

                        if (dict["$(gid)g$(nn)"].cord == 's')
                            dict["$(gid)g$(nn)"].kind = south.kind
                            dict["$(gid)g$(nn)"].cord = south.cord
                            dict["$(gid)g$(nn)"].ρ = south.ρ
                            dict["$(gid)g$(nn)"].Γ = south.Γ
                            dict["$(gid)g$(nn)"].a = south.a
                            dict["$(gid)g$(nn)"].b = south.b
                            dict["$(gid)g$(nn)"].c = south.c
                            dict["$(gid)g$(nn)"].eval = south.eval
                        end

                        if (dict["$(gid)g$(nn)"].cord == 'n')
                            dict["$(gid)g$(nn)"].kind = north.kind
                            dict["$(gid)g$(nn)"].cord = north.cord
                            dict["$(gid)g$(nn)"].ρ = north.ρ
                            dict["$(gid)g$(nn)"].Γ = north.Γ
                            dict["$(gid)g$(nn)"].a = north.a
                            dict["$(gid)g$(nn)"].b = north.b
                            dict["$(gid)g$(nn)"].c = north.c
                            dict["$(gid)g$(nn)"].eval = north.eval
                        end
                    end
                end
            end
        end

    elseif !threads
        for i in 1:mesh.l1
            for j in 1:mesh.m1
                if phi.bounds[i,j]
                    for nn in 1:phi.nbounds[i,j]
                        gid = phi.gIndex[i,j]
                        if (dict["$(gid)g$(nn)"].cord == 'w')
                            dict["$(gid)g$(nn)"].kind = west.kind
                            dict["$(gid)g$(nn)"].cord = west.cord
                            dict["$(gid)g$(nn)"].ρ = west.ρ
                            dict["$(gid)g$(nn)"].Γ = west.Γ
                            dict["$(gid)g$(nn)"].a = west.a
                            dict["$(gid)g$(nn)"].b = west.b
                            dict["$(gid)g$(nn)"].c = west.c
                            dict["$(gid)g$(nn)"].eval = west.eval
                        end

                        if (dict["$(gid)g$(nn)"].cord == 'e')
                            dict["$(gid)g$(nn)"].kind = east.kind
                            dict["$(gid)g$(nn)"].cord = east.cord
                            dict["$(gid)g$(nn)"].ρ = east.ρ
                            dict["$(gid)g$(nn)"].Γ = east.Γ
                            dict["$(gid)g$(nn)"].a = east.a
                            dict["$(gid)g$(nn)"].b = east.b
                            dict["$(gid)g$(nn)"].c = east.c
                            dict["$(gid)g$(nn)"].eval = east.eval
                        end

                        if (dict["$(gid)g$(nn)"].cord == 's')
                            dict["$(gid)g$(nn)"].kind = south.kind
                            dict["$(gid)g$(nn)"].cord = south.cord
                            dict["$(gid)g$(nn)"].ρ = south.ρ
                            dict["$(gid)g$(nn)"].Γ = south.Γ
                            dict["$(gid)g$(nn)"].a = south.a
                            dict["$(gid)g$(nn)"].b = south.b
                            dict["$(gid)g$(nn)"].c = south.c
                            dict["$(gid)g$(nn)"].eval = south.eval
                        end

                        if (dict["$(gid)g$(nn)"].cord == 'n')
                            dict["$(gid)g$(nn)"].kind = north.kind
                            dict["$(gid)g$(nn)"].cord = north.cord
                            dict["$(gid)g$(nn)"].ρ = north.ρ
                            dict["$(gid)g$(nn)"].Γ = north.Γ
                            dict["$(gid)g$(nn)"].a = north.a
                            dict["$(gid)g$(nn)"].b = north.b
                            dict["$(gid)g$(nn)"].c = north.c
                            dict["$(gid)g$(nn)"].eval = north.eval
                        end
                    end
                end
            end
        end

    end

    return dict
end

function assign_bounds!(
    phi::CSPhi3D,
    mesh::UnionCSMesh3D,
    west::BoundsStructured,
    east::BoundsStructured,
    south::BoundsStructured,
    north::BoundsStructured,
    bot::BoundsStructured,
    top::BoundsStructured,
    dict::Dict{String,BoundsStructured};
    T::Type{<:AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
    threads = false,
)
    if threads
        Base.Threads.@threads for i in 1:mesh.l1
            for j in 1:mesh.m1
                for k in 1:mesh.n1
                    if phi.bounds[i,j,k]
                        for nn in 1:phi.nbounds[i,j,k]
                            gid = phi.gIndex[i,j,k]
                            if (dict["$(gid)g$(nn)"].cord == 'w')
                                dict["$(gid)g$(nn)"].kind = west.kind
                                dict["$(gid)g$(nn)"].cord = west.cord
                                dict["$(gid)g$(nn)"].ρ = west.ρ
                                dict["$(gid)g$(nn)"].Γ = west.Γ
                                dict["$(gid)g$(nn)"].a = west.a
                                dict["$(gid)g$(nn)"].b = west.b
                                dict["$(gid)g$(nn)"].c = west.c
                                dict["$(gid)g$(nn)"].eval = west.eval
                            end

                            if (dict["$(gid)g$(nn)"].cord == 'e')
                                dict["$(gid)g$(nn)"].kind = east.kind
                                dict["$(gid)g$(nn)"].cord = east.cord
                                dict["$(gid)g$(nn)"].ρ = east.ρ
                                dict["$(gid)g$(nn)"].Γ = east.Γ
                                dict["$(gid)g$(nn)"].a = east.a
                                dict["$(gid)g$(nn)"].b = east.b
                                dict["$(gid)g$(nn)"].c = east.c
                                dict["$(gid)g$(nn)"].eval = east.eval
                            end

                            if (dict["$(gid)g$(nn)"].cord == 's')
                                dict["$(gid)g$(nn)"].kind = south.kind
                                dict["$(gid)g$(nn)"].cord = south.cord
                                dict["$(gid)g$(nn)"].ρ = south.ρ
                                dict["$(gid)g$(nn)"].Γ = south.Γ
                                dict["$(gid)g$(nn)"].a = south.a
                                dict["$(gid)g$(nn)"].b = south.b
                                dict["$(gid)g$(nn)"].c = south.c
                                dict["$(gid)g$(nn)"].eval = south.eval
                            end

                            if (dict["$(gid)g$(nn)"].cord == 'n')
                                dict["$(gid)g$(nn)"].kind = north.kind
                                dict["$(gid)g$(nn)"].cord = north.cord
                                dict["$(gid)g$(nn)"].ρ = north.ρ
                                dict["$(gid)g$(nn)"].Γ = north.Γ
                                dict["$(gid)g$(nn)"].a = north.a
                                dict["$(gid)g$(nn)"].b = north.b
                                dict["$(gid)g$(nn)"].c = north.c
                                dict["$(gid)g$(nn)"].eval = north.eval
                            end

                            if (dict["$(gid)g$(nn)"].cord == 'b')
                                dict["$(gid)g$(nn)"].kind = bot.kind
                                dict["$(gid)g$(nn)"].cord = bot.cord
                                dict["$(gid)g$(nn)"].ρ = bot.ρ
                                dict["$(gid)g$(nn)"].Γ = bot.Γ
                                dict["$(gid)g$(nn)"].a = bot.a
                                dict["$(gid)g$(nn)"].b = bot.b
                                dict["$(gid)g$(nn)"].c = bot.c
                                dict["$(gid)g$(nn)"].eval = bot.eval
                            end

                            if (dict["$(gid)g$(nn)"].cord == 't')
                                dict["$(gid)g$(nn)"].kind = top.kind
                                dict["$(gid)g$(nn)"].cord = top.cord
                                dict["$(gid)g$(nn)"].ρ = top.ρ
                                dict["$(gid)g$(nn)"].Γ = top.Γ
                                dict["$(gid)g$(nn)"].a = top.a
                                dict["$(gid)g$(nn)"].b = top.b
                                dict["$(gid)g$(nn)"].c = top.c
                                dict["$(gid)g$(nn)"].eval = top.eval
                            end
                        end
                    end
                end
            end
        end

    elseif !threads
        Base.Threads.@threads for i in 1:mesh.l1
            for j in 1:mesh.m1
                for k in 1:mesh.n1
                    if phi.bounds[i,j,k]
                        for nn in 1:phi.nbounds[i,j,k]
                            gid = phi.gIndex[i,j,k]
                            if (dict["$(gid)g$(nn)"].cord == 'w')
                                dict["$(gid)g$(nn)"].kind = west.kind
                                dict["$(gid)g$(nn)"].cord = west.cord
                                dict["$(gid)g$(nn)"].ρ = west.ρ
                                dict["$(gid)g$(nn)"].Γ = west.Γ
                                dict["$(gid)g$(nn)"].a = west.a
                                dict["$(gid)g$(nn)"].b = west.b
                                dict["$(gid)g$(nn)"].c = west.c
                                dict["$(gid)g$(nn)"].eval = west.eval
                            end

                            if (dict["$(gid)g$(nn)"].cord == 'e')
                                dict["$(gid)g$(nn)"].kind = east.kind
                                dict["$(gid)g$(nn)"].cord = east.cord
                                dict["$(gid)g$(nn)"].ρ = east.ρ
                                dict["$(gid)g$(nn)"].Γ = east.Γ
                                dict["$(gid)g$(nn)"].a = east.a
                                dict["$(gid)g$(nn)"].b = east.b
                                dict["$(gid)g$(nn)"].c = east.c
                                dict["$(gid)g$(nn)"].eval = east.eval
                            end

                            if (dict["$(gid)g$(nn)"].cord == 's')
                                dict["$(gid)g$(nn)"].kind = south.kind
                                dict["$(gid)g$(nn)"].cord = south.cord
                                dict["$(gid)g$(nn)"].ρ = south.ρ
                                dict["$(gid)g$(nn)"].Γ = south.Γ
                                dict["$(gid)g$(nn)"].a = south.a
                                dict["$(gid)g$(nn)"].b = south.b
                                dict["$(gid)g$(nn)"].c = south.c
                                dict["$(gid)g$(nn)"].eval = south.eval
                            end

                            if (dict["$(gid)g$(nn)"].cord == 'n')
                                dict["$(gid)g$(nn)"].kind = north.kind
                                dict["$(gid)g$(nn)"].cord = north.cord
                                dict["$(gid)g$(nn)"].ρ = north.ρ
                                dict["$(gid)g$(nn)"].Γ = north.Γ
                                dict["$(gid)g$(nn)"].a = north.a
                                dict["$(gid)g$(nn)"].b = north.b
                                dict["$(gid)g$(nn)"].c = north.c
                                dict["$(gid)g$(nn)"].eval = north.eval
                            end

                            if (dict["$(gid)g$(nn)"].cord == 'b')
                                dict["$(gid)g$(nn)"].kind = bot.kind
                                dict["$(gid)g$(nn)"].cord = bot.cord
                                dict["$(gid)g$(nn)"].ρ = bot.ρ
                                dict["$(gid)g$(nn)"].Γ = bot.Γ
                                dict["$(gid)g$(nn)"].a = bot.a
                                dict["$(gid)g$(nn)"].b = bot.b
                                dict["$(gid)g$(nn)"].c = bot.c
                                dict["$(gid)g$(nn)"].eval = bot.eval
                            end

                            if (dict["$(gid)g$(nn)"].cord == 't')
                                dict["$(gid)g$(nn)"].kind = top.kind
                                dict["$(gid)g$(nn)"].cord = top.cord
                                dict["$(gid)g$(nn)"].ρ = top.ρ
                                dict["$(gid)g$(nn)"].Γ = top.Γ
                                dict["$(gid)g$(nn)"].a = top.a
                                dict["$(gid)g$(nn)"].b = top.b
                                dict["$(gid)g$(nn)"].c = top.c
                                dict["$(gid)g$(nn)"].eval = top.eval
                            end
                        end
                    end
                end
            end
        end

    end

    return nothing
end
