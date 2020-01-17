# <img src="docs/src/assets/logo.svg" width="30%" align="right" /> GenomicFeatures

[![latest release][release-img]][release-url]
[![MIT license][license-img]][license-url]
[![stable documentation][docs-stable-img]][docs-stable-url]
[![latest documentation][docs-dev-img]][docs-dev-url]
![lifecycle][lifecycle-maturing]
[![Chat on Gitter][gitter-img]][gitter-url]

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

## Testing

GenomicFeatures is tested against Julia `1.X` on Linux, OS X, and Windows.

**Latest build status:**

[![travis][travis-img]][travis-url]
[![appveyor][appveyor-img]][appveyor-url]
[![coverage][codecov-img]][codecov-url]


## Contributing

We appreciate contributions from users including reporting bugs, fixing issues, improving performance and adding new features.

Take a look at the [contributing files](https://github.com/BioJulia/Contributing) detailed contributor and maintainer guidelines, and code of conduct.


### Financial contributions

We also welcome financial contributions in full transparency on our [open collective](https://opencollective.com/biojulia).
Anyone can file an expense. If the expense makes sense for the development of the community, it will be "merged" in the ledger of our open collective by the core contributors and the person who filed the expense will be reimbursed.


## Backers & Sponsors

Thank you to all our backers and sponsors!

Love our work and community? [Become a backer](https://opencollective.com/biojulia#backer).

[![backers](https://opencollective.com/biojulia/backers.svg?width=890)](https://opencollective.com/biojulia#backers)

Does your company use BioJulia?
Help keep BioJulia feature rich and healthy by [sponsoring the project](https://opencollective.com/biojulia#sponsor).
Your logo will show up here with a link to your website.

[![](https://opencollective.com/biojulia/sponsor/0/avatar.svg)](https://opencollective.com/biojulia/sponsor/0/website)
[![](https://opencollective.com/biojulia/sponsor/1/avatar.svg)](https://opencollective.com/biojulia/sponsor/1/website)
[![](https://opencollective.com/biojulia/sponsor/2/avatar.svg)](https://opencollective.com/biojulia/sponsor/2/website)
[![](https://opencollective.com/biojulia/sponsor/3/avatar.svg)](https://opencollective.com/biojulia/sponsor/3/website)
[![](https://opencollective.com/biojulia/sponsor/4/avatar.svg)](https://opencollective.com/biojulia/sponsor/4/website)
[![](https://opencollective.com/biojulia/sponsor/5/avatar.svg)](https://opencollective.com/biojulia/sponsor/5/website)
[![](https://opencollective.com/biojulia/sponsor/6/avatar.svg)](https://opencollective.com/biojulia/sponsor/6/website)
[![](https://opencollective.com/biojulia/sponsor/7/avatar.svg)](https://opencollective.com/biojulia/sponsor/7/website)
[![](https://opencollective.com/biojulia/sponsor/8/avatar.svg)](https://opencollective.com/biojulia/sponsor/8/website)
[![](https://opencollective.com/biojulia/sponsor/9/avatar.svg)](https://opencollective.com/biojulia/sponsor/9/website)


## Questions?

If you have a question about contributing or using BioJulia software, come
on over and chat to us on [Gitter](https://gitter.im/BioJulia/General), or you can try the
[Bio category of the Julia discourse site](https://discourse.julialang.org/c/domain/bio).


[release-img]:            https://img.shields.io/github/release/BioJulia/GenomicFeatures.jl.svg
[release-url]:            https://github.com/BioJulia/GenomicFeatures.jl/releases/latest
[license-img]:            https://img.shields.io/badge/license-MIT-green.svg
[license-url]:            https://github.com/BioJulia/GenomicFeatures.jl/blob/master/LICENSE
[docs-stable-img]:        https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]:        https://biojulia.github.io/GenomicFeatures.jl/stable
[docs-dev-img]:           https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]:           https://biojulia.github.io/GenomicFeatures.jl/dev/
[lifecycle-experimental]: https://img.shields.io/badge/lifecycle-experimental-orange.svg
[lifecycle-maturing]:     https://img.shields.io/badge/lifecycle-maturing-blue.svg
[lifecycle-stable]:       https://img.shields.io/badge/lifecycle-stable-brightgreen.svg
[lifecycle-retired]:      https://img.shields.io/badge/lifecycle-retired-orange.svg
[lifecycle-archived]:     https://img.shields.io/badge/lifecycle-archived-red.svg
[lifecycle-dormant]:      https://img.shields.io/badge/lifecycle-dormant-blue.svg
[lifecycle-questioning]:  https://img.shields.io/badge/lifecycle-questioning-blue.svg
[gitter-img]:             https://img.shields.io/gitter/room/BioJulia/GenomicFeatures.svg
[gitter-url]:             https://gitter.im/BioJulia/GenomicFeatures.jl
[juliapkg06-img]:         http://pkg.julialang.org/badges/GenomicFeatures_0.6.svg
[juliapkg07-img]:         http://pkg.julialang.org/badges/GenomicFeatures_0.7.svg
[juliapkg-url]:           http://pkg.julialang.org/?pkg=GenomicFeatures
[travis-img]:             https://img.shields.io/travis/BioJulia/GenomicFeatures.jl/master.svg?label=Linux+/+macOS
[travis-url]:             https://travis-ci.org/BioJulia/GenomicFeatures.jl
[appveyor-img]:           https://ci.appveyor.com/api/projects/status/dnup6vbbvai92bl8/branch/master?svg=true
[appveyor-url]:           https://ci.appveyor.com/project/BenJWard/genomicfeatures-jl/branch/master
[codecov-img]:            http://codecov.io/github/BioJulia/GenomicFeatures.jl/coverage.svg?branch=master
[codecov-url]:            http://codecov.io/github/BioJulia/GenomicFeatures.jl?branch=master
