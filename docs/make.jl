using Documenter, GenomicFeatures

DocMeta.setdocmeta!(GenomicFeatures, :DocTestSetup, :(using GenomicFeatures); recursive = true)

makedocs(
    modules = [GenomicFeatures],
    repo = Documenter.Remotes.GitHub("BioJulia", "GenomicFeatures.jl"),
    sitename = "GenomicFeatures.jl",
    doctest = true,
    linkcheck = true,
    # A workflow's page exists once the workflow is on GitHub, which a pull request that adds it precedes.
    linkcheck_ignore = [r"^https://github\.com/BioJulia/GenomicFeatures\.jl/actions/workflows/"],
    pages = [
        "Home" => "index.md",
        "Intervals" => "man/intervals.md",
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
