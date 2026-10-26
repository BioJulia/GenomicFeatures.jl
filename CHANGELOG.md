# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- CI, after BioJuliaTemplate.jl: unit tests on Julia 1 on Linux, macOS and Windows and on Julia 1.10 on Linux, with coverage, documentation with Documenter 1, Runic formatting, and Dependabot in place of CompatHelper, with the actions pinned to commits; the tests and the documentation are projects of the package's workspace.
- Gitflow, with the checks of BioJulia/Workflows: checks that pull requests into `master` come from `release/*` branches of `develop` or `hotfix/*` branches of `master`, take the next version, describe it in the changelog with its comparison links, and, if they break the public API or `master`'s tests, or narrow a compat bound, take a breaking version and say so in the changelog; and that pull requests into `develop` don't come from `release/*` or `hotfix/*` branches, and add a changelog entry.
- Back-merge: once a release or hotfix merges into `master`, the workflow of BioJulia/Workflows opens the pull request that brings `master` back into `develop`, from a `back-merge/vX.Y.Z` branch.
- This changelog, again, with every release.
- Function names in snake case, after the Julia style guide: the interval interface, `group_name`, `has_group_name`, `left_position`, `has_left_position`, `right_position`, `has_right_position`, `is_overlapping` and `is_filled`, which are BioGenerics' functions under their snake case names, so that a method defined through either name reaches callers of both; and GenomicFeatures' own `each_overlap`, `has_intersection`, `interval_type`, `base_interval_type` and `metadata_type`.

### Changed

- `coverage!` writes into a collection of intervals of any integer metadata type, `GenomicIntervalCollection{GenomicInterval{N}}` with `N <: Integer`, rather than only `UInt32`.

### Deprecated

- The squashed names of 3.0, such as `groupname`, `leftposition`, `eachoverlap` and `metadatatype`, and the internal `isordered`, which warn when deprecation warnings are on and go in 4.0; `seqname` and `hasseqname` are deprecated in favour of `group_name` and `has_group_name`.

## [3.0.1] - 2026-10-26

### Fixed

- `eachoverlap` between a stream of intervals and a `GenomicIntervalCollection` ends: a query whose group the collection lacks no longer stops the stream advancing, and a query without overlaps no longer has the rest of the stream searched again ([#52](https://github.com/BioJulia/GenomicFeatures.jl/pull/52)).

## [3.0.0] - 2023-05-01

### Breaking

- `Interval` is now `GenomicInterval`, and `IntervalCollection` is now `GenomicIntervalCollection`.
- Intervals take their group name through `groupname`, rather than `seqname`.

### Added

- `GenomicPosition`, an interval of a single position.
- The abstract types `AbstractGenomicInterval` and `AbstractGenomicCollection`, with equivalence, `isless` and `show` methods for any `AbstractGenomicInterval`.
- `groupname` and `hasgroupname` methods, exported.
- `metadatatype` and `baseintervaltype` methods, and `eltype` for generators of intervals.

### Changed

- IntervalTrees compat is `1.1`.
- The doctests run in the Documentation workflow.

## [2.1.0] - 2023-04-20

### Added

- `hasintersection` methods.
- An `eachoverlap` iterator that takes its parameters in reverse order.
- `span` and `volume` methods.
- The DOI in the README.

### Changed

- `coverage` uses `leftposition` and `rightposition`.
- Looser argument types.
- Updated workflows.

## [2.0.5] - 2021-11-27

### Added

- Downstream tests of the packages that depend on GenomicFeatures.

### Changed

- Intervals are handled through their accessor functions.
- Tidier method signatures, and the `OverlapIterator` types renamed.

### Removed

- The changelog, in favour of GitHub's release notes.

## [2.0.4] - 2021-02-20

### Changed

- DataStructures compat admits `0.18` ([#37](https://github.com/BioJulia/GenomicFeatures.jl/pull/37)).

## [2.0.3] - 2020-06-13

### Added

- Julia LTS Support
- Benchmarks

### Changed

- Documentation.
- Updated CI for General Repository.

## [2.0.2] - 2020-05-21

### Fixed

- Correct behaviour of `eachoverlap` iterator for non-overlapping queries.

## [2.0.1] - 2020-05-21

### Changed

- Documentation.
- Updated CI.

## [2.0.0] - 2020-02-16

### Breaking

- Moved BED submodule to [BED.jl](https://github.com/BioJulia/BED.jl).
- Moved BigBed submodule to [BigBed.jl](https://github.com/BioJulia/BigBed.jl).
- Moved BigWig submodule to [BigWig.jl](https://github.com/BioJulia/BigWig.jl).
- Moved GFF3 submodule to [GFF3.jl](https://github.com/BioJulia/GFF3.jl).
- Moved Indexes submodule to [Indexes.jl](https://github.com/BioJulia/Indexes.jl).
- Moved BBI submodule to [BBI.jl](https://github.com/BioJulia/BBI.jl).
- Support for julia v0.7 and v1.0 was dropped.

### Added

- BioJulia registry.

### Changed

- Migrated from BioCore to [BioGenerics](https://github.com/BioJulia/BioGenerics.jl/tree/v0.1.0).
- Updated CI.

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

## [1.0.1] - 2020-01-04 [YANKED]

### Changed

- General housekeeping.

## [1.0.0] - 2018-09-10

### Breaking

- Dropped support for julia prior to v0.7.

### Added

- Issue and PR templates.
- Code of Conduct and Contributing files.
- Support for julia v0.7 and v1.0.

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

[Unreleased]: https://github.com/BioJulia/GenomicFeatures.jl/compare/v3.0.1...HEAD
[3.0.1]: https://github.com/BioJulia/GenomicFeatures.jl/compare/v3.0.0...v3.0.1
[3.0.0]: https://github.com/BioJulia/GenomicFeatures.jl/compare/v2.1.0...v3.0.0
[2.1.0]: https://github.com/BioJulia/GenomicFeatures.jl/compare/v2.0.5...v2.1.0
[2.0.5]: https://github.com/BioJulia/GenomicFeatures.jl/compare/v2.0.4...v2.0.5
[2.0.4]: https://github.com/BioJulia/GenomicFeatures.jl/compare/v2.0.3...v2.0.4
[2.0.3]: https://github.com/BioJulia/GenomicFeatures.jl/compare/v2.0.2...v2.0.3
[2.0.2]: https://github.com/BioJulia/GenomicFeatures.jl/compare/v2.0.1...v2.0.2
[2.0.1]: https://github.com/BioJulia/GenomicFeatures.jl/compare/v2.0.0...v2.0.1
[2.0.0]: https://github.com/BioJulia/GenomicFeatures.jl/compare/v1.0.4...v2.0.0
[1.0.4]: https://github.com/BioJulia/GenomicFeatures.jl/compare/v1.0.3...v1.0.4
[1.0.3]: https://github.com/BioJulia/GenomicFeatures.jl/compare/v1.0.2...v1.0.3
[1.0.2]: https://github.com/BioJulia/GenomicFeatures.jl/compare/v1.0.1...v1.0.2
[1.0.1]: https://github.com/BioJulia/GenomicFeatures.jl/compare/v1.0.0...v1.0.1
[1.0.0]: https://github.com/BioJulia/GenomicFeatures.jl/compare/v0.2.1...v1.0.0
[0.2.1]: https://github.com/BioJulia/GenomicFeatures.jl/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/BioJulia/GenomicFeatures.jl/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/BioJulia/GenomicFeatures.jl/releases/tag/v0.1.0
