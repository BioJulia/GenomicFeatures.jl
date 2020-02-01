# GenomicFeatures.jl

[![latest release](https://img.shields.io/github/release/BioJulia/GenomicFeatures.jl.svg)](https://github.com/BioJulia/GenomicFeatures.jl/releases/latest)
[![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/BioJulia/GenomicFeatures.jl/blob/master/LICENSE)
![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)
[![Chat on Gitter](https://img.shields.io/gitter/room/BioJulia/GenomicFeatures.svg)](https://gitter.im/BioJulia/GenomicFeatures.jl)

## Description
GenomicFeatures provides utilities for working with interval based genomic annotations.
It builds on [IntervalTrees](https://github.com/biojulia/intervaltrees.jl) to provide a data-structures and algorithms for various formats such as [BED](https://github.com/biojulia/bed.jl), [GFF3](https://github.com/biojulia/gff3.jl), [bigWig](https://github.com/biojulia/bigwig.jl) and [bigBed](https://github.com/biojulia/bigbed.jl).

## Installation
Releases of GenomicFeatures version 2.0.0 and above are registered and made available to install through BioJulia's package registry.
By default, Julia's package manager only uses the "General" package registry.
Your Julia configuration needs to include the BioJulia registry to be able to install the latest version of GenomicFeatures.

To add the BioJulia registry from the [Julia REPL](https://docs.julialang.org/en/v1/manual/getting-started/), press `]` to enter [pkg mode](https://docs.julialang.org/en/v1/stdlib/Pkg/), then enter the following command:
```julia
registry add https://github.com/BioJulia/BioJuliaRegistry.git
```

After adding the registry to your configuration, you can install GenomicFeatures while in [pkg mode](https://docs.julialang.org/en/v1/stdlib/Pkg/) with the following:
```julia
add GenomicFeatures
```

If you are interested in the cutting edge of the development, please check out the [develop branch](https://github.com/BioJulia/GenomicFeatures.jl/tree/develop) to try new features before release.
