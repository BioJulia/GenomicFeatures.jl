using Pkg

Pkg.activate(@__DIR__)
Pkg.instantiate()

Pkg.status()

using PkgBenchmark

results = benchmarkpkg(
    dirname(@__DIR__),
    BenchmarkConfig(
        env = Dict(
            "JULIA_NUM_THREADS" => "1",
            "OMP_NUM_THREADS" => "1",
        ),
    )
)

dir_results = joinpath(@__DIR__, "results")
mkpath(dir_results)

writeresults(joinpath(dir_results, "$(results.date).json"), results)
export_markdown(joinpath(dir_results, "$(results.date).md"), results)
