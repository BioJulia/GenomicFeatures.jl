# The regression suite, which AirspeedVelocity.jl runs at each revision to compare their times and allocations; see docs/src/benchmarks.md.
# AirspeedVelocity runs a copy of this script away from the repository, so it is self-contained.
using BenchmarkTools
using GenomicFeatures
using Random

# `n` random intervals over three groups of `span` positions each, as long as chromosomes, of 1 to 2000 positions, with their index as metadata.
function random_intervals(n; seed = 1234, groups = ["chr1", "chr2", "chr3"], span = 100_000_000)
    rng = MersenneTwister(seed)
    intervals = Vector{GenomicInterval{Int}}(undef, n)
    for i in 1:n
        len = rand(rng, 1:2000)
        first = rand(rng, 1:(span - len + 1))
        intervals[i] = GenomicInterval{Int}(rand(rng, groups), first, first + len - 1, rand(rng, (STRAND_POS, STRAND_NEG)), i)
    end
    return intervals
end

const SIZES = (1_000, 100_000)

const SUITE = BenchmarkGroup()
for name in ("accessors", "sort", "collection", "each_overlap", "coverage", "has_intersection")
    SUITE[name] = BenchmarkGroup()
end

for n in SIZES
    intervals = random_intervals(n)
    sorted = sort(intervals)
    others = sort(random_intervals(n; seed = 4321))
    collection = GenomicIntervalCollection(sorted)
    other_collection = GenomicIntervalCollection(others)
    label = string(n)

    accessors = SUITE["accessors"][label] = BenchmarkGroup()
    accessors["group_name"] = @benchmarkable group_name.($intervals)
    accessors["left_position"] = @benchmarkable left_position.($intervals)
    accessors["right_position"] = @benchmarkable right_position.($intervals)
    accessors["strand"] = @benchmarkable strand.($intervals)
    accessors["metadata"] = @benchmarkable metadata.($intervals)

    SUITE["sort"][label] = @benchmarkable sort!(unsorted) setup = (unsorted = copy($intervals)) evals = 1

    building = SUITE["collection"][label] = BenchmarkGroup()
    building["sorted vector"] = @benchmarkable GenomicIntervalCollection($sorted)
    building["push!"] = @benchmarkable foreach(interval -> push!(empty, interval), $intervals) setup = (empty = GenomicIntervalCollection{GenomicInterval{Int}}()) evals = 1

    overlaps = SUITE["each_overlap"][label] = BenchmarkGroup()
    overlaps["vector, vector"] = @benchmarkable collect(each_overlap($sorted, $others))
    overlaps["vector, collection"] = @benchmarkable collect(each_overlap($sorted, $other_collection))
    overlaps["collection, collection"] = @benchmarkable collect(each_overlap($collection, $other_collection))

    SUITE["coverage"][label] = @benchmarkable coverage($collection)

    SUITE["has_intersection"][label] = @benchmarkable count(has_intersection($collection), $others)
end
