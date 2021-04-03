"""
This module provides:
    -

The exported names are:
$(EXPORTS)
"""
module VisualizationData

using DocStringExtensions
using WriteVTK

using dopingFVM.Structures
using dopingFVM.Tools

export check_folder
export plot_paraviewVTK
export plot_paraviewPVD

include("CheckFolder.jl")

include("CSParaviewVTK.jl")

include("CSParaviewPVD.jl")

end # module
