using Documenter, GenomicFeatures

makedocs()
deploydocs(
    deps = Deps.pip("mkdocs", "pygments", "mkdocs-material"),
    repo = "github.com/BioJulia/GenomicFeatures.jl.git",
    julia = "0.5",
    osname = "linux",
)
