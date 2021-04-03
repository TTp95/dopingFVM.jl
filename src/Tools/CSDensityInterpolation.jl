"""

"""
@inline function density_interpolation(
    lenghts1::AbstractFloat,
    lenghts2::AbstractFloat,
    values1::AbstractFloat,
    values2::AbstractFloat;
    interpolation = 1,
)
    #Inicializate value
    value = 0.0

    #Interpolation mode selection
    if (interpolation == 1)
        num = values1 * (0.5 * lenghts2) + values2 * (0.5 * lenghts1)
        den = 0.5 * (lenghts1 + lenghts2)
        value = num / den
    else
        error("density_interpolation mode unimplemented")
    end

    return value
end
