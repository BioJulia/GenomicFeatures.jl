using BenchmarkTools
# using BenchmarkTools: @benchmarkable, BenchmarkGroup

using GenomicFeatures

include(joinpath(@__DIR__, "..", "test", "Utilities.jl"))

import ..Utilities: random_intervals

N = 1000
SEED = 1234
SEQNAMES = "chr" .* string.(1:3)

intervals = random_intervals(SEQNAMES, 1000, N, SEED)
intervals_sorted = sort(intervals)

SUITE = BenchmarkGroup()

let suite = SUITE["accessors"] = BenchmarkGroup()
    s0 = suite["$(typeof(intervals))"] = BenchmarkGroup()
    s0["seqname"] = @benchmarkable(seqname.($intervals))
    s0["leftposition"] = @benchmarkable(leftposition.($intervals))
    s0["rightposition"] = @benchmarkable(rightposition.($intervals))
    s0["strand"] = @benchmarkable(strand.($intervals))
    s0["metadata"] = @benchmarkable(metadata.($intervals))
end

let suite = SUITE["sort"] = BenchmarkGroup()
    suite["$(typeof(intervals))"] = @benchmarkable(sort(i), setup=(i = copy($intervals)))
end

let suite = SUITE["insert"] = BenchmarkGroup()
    suite["shorthand"] = @benchmarkable(GenomicIntervalCollection($intervals_sorted))
    suite["type"] = @benchmarkable(GenomicIntervalCollection{Int}($intervals_sorted))
end

let suite = SUITE["push"] = BenchmarkGroup()
    suite["$(typeof(intervals))"] = @benchmarkable([push!(col, i) for i in $intervals], setup=(col=GenomicIntervalCollection{Int}()))
end

let suite = SUITE["eachoverlap"] = BenchmarkGroup()
    intervals_a = intervals_sorted
    intervals_b = sort(random_intervals(SEQNAMES, 1000, N, SEED+1))

    col_a = GenomicIntervalCollection(intervals_a)
    col_b = GenomicIntervalCollection(intervals_b)

    As = [intervals_a, col_a]
    Bs = [intervals_b, col_b]

    for (A, B) in Iterators.product(As,Bs)
        str = "$(typeof(A)), $(typeof(B))"
        suite[str] = @benchmarkable(collect(eachoverlap($A,$B)))
    end
end
