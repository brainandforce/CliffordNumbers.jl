using CliffordNumbers
using Documenter

DocMeta.setdocmeta!(CliffordNumbers, :DocTestSetup, :(using CliffordNumbers); recursive=true)

is_ci_env = (get(ENV, "CI", nothing) == true)
@info "is_ci_env == $is_ci_env"

makedocs(;
    sitename="CliffordNumbers.jl",
    modules=[CliffordNumbers],
    checkdocs = :exports,
    format=Documenter.HTML(;
        prettyurls = is_ci_env,
        canonical = "https://brainandforce.github.io/CliffordNumbers.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages=[
        "Home" => "index.md",
        "Operations" => "operations.md",
        "API" => "api.md"
    ],
)

deploydocs(;
    repo="github.com/brainandforce/CliffordNumbers.jl",
    devbranch="main",
)
