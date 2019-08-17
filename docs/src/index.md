GenomicFeatures.jl
==================

## Description
GenomicFeatures provides utilities for working with interval based genomic annotations.
It builds on [IntervalTrees](https://github.com/biojulia/intervaltrees.jl) to provide a data-structures and algorithms for various formats such as [BED](https://github.com/biojulia/bed.jl), [GFF3](https://github.com/biojulia/gff3.jl), [bigWig](https://github.com/biojulia/bigwig.jl) and [bigBed](https://github.com/biojulia/bigbed.jl).

## Installation
GenomicFeatures is made available to install through BioJulia's package registry.

Julia by default only watches the "General" package registry, so before you start, you should add the BioJulia package registry.

Start a julia terminal, hit the `]` key to enter pkg mode (you should see the prompt change from `julia>` to `pkg>` ), then enter the following command:

```julia
registry add https://github.com/BioJulia/BioJuliaRegistry.git
```

After you've added the registry, you can install GenomicFeatures from the julia REPL.
Press `]` to enter pkg mode again, and enter the following:

```julia
add GenomicFeatures
```

If you are interested in the cutting edge of the development, please check out the [develop branch](https://github.com/BioJulia/GenomicFeatures.jl/tree/develop) to try new features before release.
