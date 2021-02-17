using dopingFVM
using Documenter

makedocs(;
    modules=[dopingFVM],
    authors="Felipe A. Díaz",
    repo="https://github.com/TTp95/dopingFVM.jl/blob/{commit}{path}#L{line}",
    sitename="dopingFVM.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://TTp95.github.io/dopingFVM.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/TTp95/dopingFVM.jl",
)
