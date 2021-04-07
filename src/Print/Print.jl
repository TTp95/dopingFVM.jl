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

include("PrintConsole.jl")

end # module
