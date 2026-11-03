using Documenter, GenomicFeatures

DocMeta.setdocmeta!(GenomicFeatures, :DocTestSetup, :(using GenomicFeatures); recursive = true)

makedocs(
    modules = [GenomicFeatures],
    repo = Documenter.Remotes.GitHub("BioJulia", "GenomicFeatures.jl"),
    sitename = "GenomicFeatures.jl",
    doctest = true,
    linkcheck = true,
    # The pages of this repository's own workflows and files on develop exist once they are on GitHub, which a pull request that adds them precedes.
    linkcheck_ignore = [r"^https://github\.com/BioJulia/GenomicFeatures\.jl/(actions/workflows|blob/develop)/"],
    pages = [
        "Home" => "index.md",
        "Intervals" => "man/intervals.md",
        "Benchmarks" => "benchmarks.md",
        "API Reference" => "api/api.md",
    ],
    format = Documenter.HTML(edit_link = "develop"),
    authors = "Kenta Sato, Sabrina J. Ward, Ciarán O'Mara, The BioJulia Organisation, and other contributors.",
    checkdocs = :public,
)

deploydocs(;
    repo = "github.com/BioJulia/GenomicFeatures.jl.git",
    devbranch = "develop",
    push_preview = true,
    deps = nothing,
    make = nothing,
)
