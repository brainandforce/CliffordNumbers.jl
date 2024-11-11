using CliffordNumbers
using Quaternions
using Documenter

import CliffordNumbers.BaseNumber

using_directives = :(using CliffordNumbers; using Quaternions)
DocMeta.setdocmeta!(CliffordNumbers, :DocTestSetup, using_directives; recursive=true)

is_ci_env = (get(ENV, "CI", nothing) == true)
@info "is_ci_env == $is_ci_env"

makedocs(;
    sitename="CliffordNumbers.jl",
    modules=[CliffordNumbers, Base.get_extension(CliffordNumbers, :CliffordNumbersQuaternionsExt)],
    doctest=false,
    checkdocs = :exports,
    format=Documenter.HTML(;
        prettyurls = is_ci_env,
        canonical = "https://brainandforce.github.io/CliffordNumbers.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages=[
        "Home" => "index.md",
        "Clifford number types" => "numeric.md",
        "Metric signatures" => "metrics.md",
        "Indexing" => "indexing.md",
        "Operations" => "operations.md",
        "Performance tips" => "performance.md",
        "Extensions" => "extensions.md",
        "API" => Any[
            "CliffordNumbers" => "api/clifford.md",
            "Indexing" => "api/indexing.md",
            "Math" => "api/math.md",
            "Metric signatures" => "api/metrics.md",
            "Internals" => "api/internal.md",
        ]
    ],
)

deploydocs(;
    repo="github.com/brainandforce/CliffordNumbers.jl.git",
    devbranch="main",
)
