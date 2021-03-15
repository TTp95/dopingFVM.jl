"""

"""
function density_interpolationCS(
    lenghts::Union{Array{<:AbstractFloat,1}, Array{<:AbstractFloat,2}},
    values::Union{Array{<:AbstractFloat,1}, Array{<:AbstractFloat,2}};
    interpolation = 1
)
    #Inicializate value
    value = 0.0

    #Interpolation mode selection
    if (interpolation == 1)
        num = values[1] * (0.5 * lenghts[2]) + values[2] * (0.5 * lenghts[1])
        den = 0.5 * (lenghts[1] + lenghts[2])
        value = num / den
    else
        error("density_interpolationCS mode unimplemented")
    end

    return value
end
