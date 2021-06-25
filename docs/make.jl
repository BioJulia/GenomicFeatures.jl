using Pkg

using Documenter, GenomicFeatures

DocMeta.setdocmeta!(GenomicFeatures, :DocTestSetup, :(using GenomicFeatures); recursive=true)

format = Documenter.HTML(
    edit_link = "develop"
)

makedocs(
    modules = [GenomicFeatures],
    checkdocs = :all,
    linkcheck = true,
    format = format,
    sitename = "GenomicFeatures.jl",
    authors = replace(join(Pkg.TOML.parsefile("Project.toml")["authors"], ", "), r" <.*?>" => "" ) * ", The BioJulia Organisation, and other contributors.",
    pages = [
        "Home" => "index.md",
        "Intervals" => "man/intervals.md",
        "API Reference" => "api/api.md"
    ],
)

deploydocs(
    repo = "github.com/BioJulia/GenomicFeatures.jl.git",
    devbranch = "develop",
    push_preview = true
)
