using Documenter
using VectorUtils

makedocs(
    sitename = "VectorUtils.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://mirajcs.github.io/VectorUtils",
        assets = String[],
    ),
    modules = [VectorUtils],
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "API Reference" => [
            "Core Functions" => "api/core.md",
            "Symbolic Computations" => "api/symbolic.md",
            "Numeric Computations" => "api/numeric.md",
        ],
        "Examples" => [
            "Basic Usage" => "examples/basic.md",
            "Parametric Curves" => "examples/curves.md",
            "Frenet Frames" => "examples/frenet.md",
        ],
        "Mathematical Background" => "theory.md",
    ],
    repo = "https://github.com/mirajcs/VectorUtils",
    authors = "Miraj Samarakkody",
)

deploydocs(
    repo = "github.com/mirajcs/VectorUtils.git",
    devbranch = "main",
)