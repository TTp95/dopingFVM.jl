"""
This module provides:
    -

The exported names are:
$(EXPORTS)
"""
module Bounds

using Base.Threads

using DocStringExtensions

using dopingFVM.Structures
using dopingFVM.StructuresConstructors
using dopingFVM.Tools

export check_bounds!
export create_BoundsDict
export bounds_template
export assign_bounds!
export _evaluate_bounds!_

include("CSBoundsCheck.jl")

include("CSBoundsDictionary.jl")

include("CSBoundsTemplate.jl")

include("CSAssignBounds.jl")

include("CSEvaluateBounds.jl")

end # module
