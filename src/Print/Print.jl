"""
This module provides:
    -

The exported names are:
$(EXPORTS)
"""
module Print

using DocStringExtensions
using Printf
using WriteVTK

using dopingFVM.Structures
using dopingFVM.Tools

export print_consoleError
export printf_consoleError
export print_summary

include("PrintConsoleError.jl")
include("PrintfConsoleError.jl")

include("PrintSummary.jl")

end # module
