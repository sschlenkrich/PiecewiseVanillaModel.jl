using PiecewiseVanillaModel
using Documenter

DocMeta.setdocmeta!(PiecewiseVanillaModel, :DocTestSetup, :(using PiecewiseVanillaModel); recursive=true)

makedocs(;
    modules=[PiecewiseVanillaModel],
    authors="Sebastian Schlenkrich",
    sitename="PiecewiseVanillaModel.jl",
    format=Documenter.HTML(;
        canonical="https://sschlenkrich.github.io/PiecewiseVanillaModel.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/sschlenkrich/PiecewiseVanillaModel.jl",
    devbranch="main",
)
