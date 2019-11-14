using Documenter, GenomicFeatures

makedocs(
    format = :html,
    sitename = "GenomicFeatures.jl",
    pages = [
        "Home" => "index.md",
        "Intervals" => "intervals.md",
    ],
    authors = "Kenta Sato, D. C. Jones, Ben J. Ward, The BioJulia Organisation and other contributors."
)
deploydocs(
    repo = "github.com/BioJulia/GenomicFeatures.jl.git",
    julia = "0.6",
    osname = "linux",
    target = "build",
    deps = nothing,
    make = nothing
)
