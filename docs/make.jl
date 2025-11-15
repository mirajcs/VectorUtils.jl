using Documenter
using VectorUtils

makedocs(
    sitename = "VectorUtils.jl",
    modules = [VectorUtils],
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md"
    ],
)