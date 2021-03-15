"""
mutable struct SystemTime{T <: AbstractFloat} <: TypeDopingFMV

"""
function SystemTime(; T::Type{<:AbstractFloat} = Float64)
    return DeltaTime(T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0))
end
