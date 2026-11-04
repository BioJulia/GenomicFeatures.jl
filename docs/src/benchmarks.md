# Benchmarks

GenomicFeatures' benchmarks are in `benchmark/`, a project of the package's workspace.

`benchmark/benchmarks.jl` defines `SUITE`, a [BenchmarkTools.jl](https://github.com/JuliaCI/BenchmarkTools.jl) group of the accessors, `sort`, building a `GenomicIntervalCollection`, `eachoverlap`, `coverage` and `hasintersection`, on 10³ and 10⁵ random intervals.
[AirspeedVelocity.jl](https://github.com/MilesCranmer/AirspeedVelocity.jl) runs it at two revisions and compares their times and allocations.
On every pull request, the Benchmarks workflow compares the pull request with the commit of its base branch, `develop` or `master`, that it merges onto, and comments the table, with the ratio of their times, on the pull request.
It runs the pull request's own suite at both revisions, so the suite uses only names that the base has too.
A maintainer's `skip benchmark` label skips them, for a change that can't affect them, such as to the documentation or the CI, and removing the label runs them; Dependabot's updates of the actions carry it.

To compare revisions locally, install AirspeedVelocity's `benchpkg` and run it from the repository, such as for `develop` and the working tree:

```sh
julia -e 'using Pkg; Pkg.add("AirspeedVelocity"); Pkg.build("AirspeedVelocity")'
~/.julia/bin/benchpkg GenomicFeatures --rev=develop,dirty --bench-on=dirty
```
