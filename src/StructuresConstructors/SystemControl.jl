"""


"""
function create_SystemControl(;
    T::Type{AbstractFloat} = Float64,
    N::Type{<:Signed} = Int64,
)
    return SystemControl{T,N}(
        "",
        "results",
        true,
        0,
        0,
        0,
        0,
        0,
        false,
        1.0e20,
        1,
        0.0,
        0.0,
        false,
        false

    )
end
