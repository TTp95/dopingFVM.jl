"""

"""
@inline function gamma_interpolation(
    lenghts1::AbstractFloat,
    lenghts2::AbstractFloat,
    values1::AbstractFloat,
    values2::AbstractFloat;
    interpolation = 2
)
    #Inicializate value
    value = 0.0

    #Interpolation mode selection
    if (interpolation == 1) #Simple interpolation
        num = values1 * (0.5 * lenghts2) + values2 * (0.5 * lenghts1)
        den = 0.5 * (lenghts1 + lenghts2)
        value = num / den

    elseif (interpolation == 2) #Interpolation heat
        g = (0.5 * lenghts1) / (0.5 * (lenghts1 + lenghts2))
        num = values2 * values1
        den = values1 * (1.0 - g) + values2 * g
        value = num / den

    else
        error("gamma_interpolation mode unimplemented")

    end

    return value
end
