"""
crete_DeltaTime(; T::Type{<:AbstractFloat} = Float64)

"""
function crete_DeltaTime(; T::Type{<:AbstractFloat} = Float64)
    return DeltaTime(T(0.0), T(0.0), T(0.0))
end
