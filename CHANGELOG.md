# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Added
- BioJulia registry.

## [2.0.0] - 2020-01-

### Changed
- Migrated from BioCore to [BioGenerics](https://github.com/BioJulia/BioGenerics.jl/tree/v0.1.0).
- Updated CI.

### Removed
- :exclamation: Moved BED submodule to [BED.jl](https://github.com/BioJulia/BED.jl).
- :exclamation: Moved BigBed submodule to [BigBed.jl](https://github.com/BioJulia/BigBed.jl).
- :exclamation: Moved BigWig submodule to [BigWig.jl](https://github.com/BioJulia/BigWig.jl).
- :exclamation: Moved GFF3 submodule to [GFF3.jl](https://github.com/BioJulia/GFF3.jl).
- :exclamation: Moved Indexes submodule to [Indexes.jl](https://github.com/BioJulia/Indexes.jl).
- :exclamation: Moved BBI submodule to [BBI.jl](https://github.com/BioJulia/BBI.jl).
- :exclamation: Support for julia v0.7 and v1.0 was dropped.

## [1.0.4] - 2020-01-15
### Changed
- General housekeeping (code).

## [1.0.3] - 2020-01-07
### Changed
- Add ColorTypes "0.9" to compat.
- Add FixedPointNumbers "0.7" to compat.

### Fixed
- CompatHelper configuration.

## [1.0.2] - 2020-01-07
### Fixed
- Corrected typo in Strand constructor.

## [1.0.1] - 2020-01-04 (yanked)
### Changed
- General housekeeping.

## [1.0.0] - 2018-09-10
### Added
- Issue and PR templates.
- Code of Conduct and Contributing files.
- Support for julia v0.7 and v1.0.

### Removed
- :exclamation: Support for julia prior to v0.7.

## [0.2.1] - 2017-06-27
### Added
- Read and write methods for Strand type ([#4](https://github.com/BioJulia/GenomicFeatures.jl/pull/4)).

## [0.2.0] - 2017-06-25
### Added
- Optional filter functions to `eachoverlap`, and a `findfirst` function ([#3](https://github.com/BioJulia/GenomicFeatures.jl/pull/3)).

## [0.1.0] - 2017-06-17
### Added
- Code from Bio.jl.
- Support for GFF3 and BigWig.

[Unreleased]: https://github.com/BioJulia/GenomicFeatures.jl/compare/v2.0.0...HEAD
[1.0.4]: https://github.com/BioJulia/GenomicFeatures.jl/compare/v1.0.4...v2.0.0
[1.0.4]: https://github.com/BioJulia/GenomicFeatures.jl/compare/v1.0.3...v1.0.4
[1.0.3]: https://github.com/BioJulia/GenomicFeatures.jl/compare/v1.0.2...v1.0.3
[1.0.2]: https://github.com/BioJulia/GenomicFeatures.jl/compare/v1.0.1...v1.0.2
[1.0.1]: https://github.com/BioJulia/GenomicFeatures.jl/compare/v1.0.0...v1.0.1
[1.0.0]: https://github.com/BioJulia/GenomicFeatures.jl/compare/v0.2.1...v1.0.0
[0.2.1]: https://github.com/BioJulia/GenomicFeatures.jl/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/BioJulia/GenomicFeatures.jl/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/BioJulia/GenomicFeatures.jl/tree/v0.1.0
