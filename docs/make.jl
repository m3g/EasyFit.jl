using Documenter
using EasyFit

makedocs(;
    modules=[EasyFit],
    authors="Leandro Martínez <lmartine@unicamp.br> and contributors",
    sitename="EasyFit.jl",
    doctest=false,
    checkdocs=:none,
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://m3g.github.io/EasyFit.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Linear fit" => "linear.md",
        "Quadratic fit" => "quadratic.md",
        "Cubic fit" => "cubic.md",
        "N-th degree polynomial fit" => "ndegree.md",
        "Exponential fits" => "exponential.md",
        "Normalized exponential decay" => "expdecay.md",
        "Splines" => "splines.md",
        "Moving averages" => "movingaverage.md",
        "Density function" => "density.md",
        "Bounds and fixed parameters" => "bounds.md",
        "Options" => "options.md",
    ],
)

deploydocs(;
    repo="github.com/m3g/EasyFit.jl",
    devbranch="main",
)
