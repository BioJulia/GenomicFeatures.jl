using Pkg

Pkg.develop(PackageSpec(path=pwd()))
Pkg.instantiate()

using Documenter, GenomicFeatures

makedocs(
    format = Documenter.HTML(
        edit_link = :commit
    ),
    modules = [GenomicFeatures],
    sitename = "GenomicFeatures.jl",
    pages = [
        "Home" => "index.md",
        "Intervals" => "man/intervals.md",
        "I/O" => [
            "BED" => "io/bed.md",
            "GFF3" => "io/gff3.md",
            "BigWig" => "io/bigwig.md",
            "BigBed" => "io/bigbed.md"
        ],
        "API Reference" => "api/api.md"
    ],
    authors = replace(join(Pkg.TOML.parsefile("Project.toml")["authors"], ", "), r" <.*?>" => "" ) * ", The BioJulia Organisation, and other contributors."
)
deploydocs(
    repo = "github.com/BioJulia/GenomicFeatures.jl.git",
    devbranch = "develop"
)
