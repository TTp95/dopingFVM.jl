"""

"""
function gamma_interpolationCS(
    lenghts::Union{Array{<:AbstractFloat,1}, Array{<:AbstractFloat,2}},
    values::Union{Array{<:AbstractFloat,1}, Array{<:AbstractFloat,2}};
    interpolation = 1
)
    #Inicializate value
    value = 0.0

    #Interpolation mode selection
    if (interpolation == 1) #Simple interpolation
        num = values[1] * (0.5 * lenghts[2]) + values[2] * (0.5 * lenghts[1])
        den = 0.5 * (lenghts[1] + lenghts[2])
        value = num / den

    elseif (interpolation == 2) #Interpolation heat
        g = (0.5 * lenghts[1]) / (0.5 * (lenghts[1] + lenghts[2]))
        num = values[2] * values[1]
        den = values[1] * (1.0 - g) + values[2] * g
        value = num / den

    else
        error("gamma_interpolationCS mode unimplemented")

    end

    return value
end
