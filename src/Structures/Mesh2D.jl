"""
    mysearch(array::MyArray{T}, val::T; verbose=true) where {T} -> Int
    mysearch(array::MyArray{T}, val::T; verbose=true) where {T} -> Int

Searches the `array` for the `val`. For some reason we don't want to use Julia's
builtin search :)

# Arguments
- `array::MyArray{T}`: the array to search
- `val::T`: the value to search for

# Keywords
- `verbose::Bool=true`: print out progress details

# Returns
- `Int`: the index where `val` is located in the `array`

# Throws
- `NotFoundError`: I guess we could throw an error if `val` isn't found.
"""

struct M
    x::Array{Float64,1}         #Posición nodo, eje x
    y::Array{Float64,1}         #Posición nodo, eje y
    dx::Array{Float64,1}        #Longitud del vc eje x
    dy::Array{Float64,1}        #Longitud del vc eje y
    vol::Array{Float64,2}       #Volumen de vc
    l1::Int64                   #Último nodo en eje x
    m1::Int64                   #Último nodo en eje y
end
