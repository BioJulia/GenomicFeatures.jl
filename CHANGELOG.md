# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Support for julia v1.1 and v1.2.
- Autodocs for public functions.
- `GenomicInterval` is now a subtype of `AbstractGenomicInterval`.

### Changed
- Migrated from BioCore to [BioGenerics](https://github.com/BioJulia/BioGenerics.jl/tree/v0.1.0).
- Renamed `Interval` to `GenomicInterval`.
- Renamed `IntervalCollection` to `GenomicIntervalCollection`.


### Removed
- :exclamation: BED module was moved to [BED.jl](https://github.com/BioJulia/BED.jl).
- :exclamation: BigBed module was moved to [BigBed.jl](https://github.com/BioJulia/BigBed.jl).
- :exclamation: BigWig module was moved to [BigWig.jl](https://github.com/BioJulia/BigWig.jl).
- :exclamation: GFF3 module was moved to [GFF3.jl](https://github.com/BioJulia/GFF3.jl).
- :exclamation: Indexes module was moved to [Indexes.jl](https://github.com/BioJulia/Indexes.jl).
- :exclamation: BBI module was moved to [BBI.jl](https://github.com/BioJulia/BBI.jl).
- :exclamation: Support for julia v0.7 and v1.0 has been dropped.

## [1.0.0]
- Version 1.0 updates ([#15](https://github.com/BioJulia/GenomicFeatures.jl/pull/15))

## [0.2.1]
### Added
- Read and write methods for Strand ([#4](https://github.com/BioJulia/GenomicFeatures.jl/pull/4)).

## [0.2.0]
### Added
- Optional filter functions to `eachoverlap`, and a `findfirst` function ([#3](https://github.com/BioJulia/GenomicFeatures.jl/pull/3)).

## [0.1.0]
- Move code from Bio.jl.
- Add support for GFF3 and BigWig.

[Unreleased]: https://github.com/BioJulia/GenomicFeatures.jl/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/BioJulia/GenomicFeatures.jl/compare/v0.2.1...v1.0.0
[0.2.1]: https://github.com/BioJulia/GenomicFeatures.jl/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/BioJulia/GenomicFeatures.jl/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/BioJulia/GenomicFeatures.jl/tree/v0.1.0
