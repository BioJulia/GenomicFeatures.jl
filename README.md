# <img src="./sticker.svg" width="30%" align="right" /> GenomicFeatures

[![latest release](https://img.shields.io/github/release/BioJulia/GenomicFeatures.jl.svg)](https://github.com/BioJulia/GenomicFeatures.jl/releases/latest)
[![DOI](https://zenodo.org/badge/94160625.svg)](https://zenodo.org/badge/latestdoi/94160625)
[![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/BioJulia/GenomicFeatures.jl/blob/master/LICENSE)
[![Stable documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://biojulia.github.io/GenomicFeatures.jl/stable)
[![Latest documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://biojulia.github.io/GenomicFeatures.jl/dev/)
[![lifecycle](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

## Description
The GenomicFeatures package provides utilities for working with interval based genomic annotations.
It builds on [IntervalTrees](https://github.com/biojulia/intervaltrees.jl) to provide a data-structures and algorithms for various formats such as [BED](https://github.com/biojulia/bed.jl), [GFF3](https://github.com/biojulia/gff3.jl), [bigWig](https://github.com/biojulia/bigwig.jl) and [bigBed](https://github.com/biojulia/bigbed.jl).  

## Installation
You can install the GenomicFeatures package from the [Julia REPL](https://docs.julialang.org/en/v1/manual/getting-started/).
Press `]` to enter [pkg mode](https://docs.julialang.org/en/v1/stdlib/Pkg/), then enter the following command:
```julia
add GenomicFeatures
```

If you are interested in the cutting edge of the development, please check out the [develop branch](https://github.com/BioJulia/GenomicFeatures.jl/tree/develop) to try new features before release.

## Testing

GenomicFeatures is tested against Julia 1 on Linux, macOS and Windows, and against Julia 1.10, the long-term support release, on Linux.

**Latest build status:**

[![Unit tests](https://github.com/BioJulia/GenomicFeatures.jl/actions/workflows/UnitTesting.yml/badge.svg)](https://github.com/BioJulia/GenomicFeatures.jl/actions/workflows/UnitTesting.yml)
[![Runic](https://github.com/BioJulia/GenomicFeatures.jl/actions/workflows/Formatting.yml/badge.svg)](https://github.com/BioJulia/GenomicFeatures.jl/actions/workflows/Formatting.yml)
[![Documentation](https://github.com/BioJulia/GenomicFeatures.jl/actions/workflows/Documentation.yml/badge.svg)](https://github.com/BioJulia/GenomicFeatures.jl/actions/workflows/Documentation.yml)
[![codecov](https://codecov.io/gh/BioJulia/GenomicFeatures.jl/branch/develop/graph/badge.svg)](https://codecov.io/gh/BioJulia/GenomicFeatures.jl)
[![Downstream](https://github.com/BioJulia/GenomicFeatures.jl/actions/workflows/Downstream.yml/badge.svg)](https://github.com/BioJulia/GenomicFeatures.jl/actions/workflows/Downstream.yml)

## Contributing

We appreciate [contributions](https://github.com/BioJulia/GenomicFeatures.jl/graphs/contributors) from users including reporting bugs, fixing issues, improving performance and adding new features.
See the [contributing files](https://github.com/BioJulia/Contributing) for detailed contributor and maintainer guidelines and the code of conduct.

## Questions

If you have a question about contributing or using BioJulia software, come and chat with us in the `#biology` channel of [the Julia Slack workspace](https://julialang.slack.com/archives/biology), or try the [Bio category of the Julia Discourse site](https://discourse.julialang.org/c/domain/bio/15).
