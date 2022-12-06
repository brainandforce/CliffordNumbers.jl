using CliffordNumbers
using Documenter

DocMeta.setdocmeta!(CliffordNumbers, :DocTestSetup, :(using CliffordNumbers); recursive=true)

makedocs(;
    modules=[CliffordNumbers],
    authors="Brandon Flores",
    repo="https://github.com/brainandforce/CliffordNumbers.jl/blob/{commit}{path}#{line}",
    sitename="CliffordNumbers.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://brainandforce.github.io/CliffordNumbers.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/brainandforce/CliffordNumbers.jl",
    devbranch="main",
)
