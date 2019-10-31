using GenomicFeatures
using Test

using Distributions
import Random


# Test that an array of intervals is well ordered
function is_all_ordered(intervals::Vector{I}) where I <: GenomicInterval
    for i = 2:length(intervals)
        if !GenomicFeatures.isordered(intervals[i-1], intervals[i])
            return false
        end
    end
    return true
end

# Generate an array of n random GenomicInterval{Int} object.
# With sequence names, samples from seqnames, and intervals drawn to lie in [1, maxpos].
function random_intervals(seqnames, maxpos::Int, n::Int)
    seq_dist = Categorical(length(seqnames))
    strand_dist = Categorical(2)
    length_dist = Normal(1000, 1000)
    intervals = Vector{GenomicInterval{Int}}(undef, n)
    for i in 1:n
        intlen = maxpos
        while intlen >= maxpos || intlen <= 0
            intlen = ceil(Int, rand(length_dist))
        end
        first = rand(1:maxpos-intlen)
        last = first + intlen - 1
        strand = rand(strand_dist) == 1 ? STRAND_POS : STRAND_NEG
        intervals[i] = GenomicInterval{Int}(seqnames[rand(seq_dist)], first, last, strand, i)
    end
    return intervals
end

# A simple interval intersection implementation to test against.
function simple_intersection(intervals_a, intervals_b; filter=(a,b)->true)
    sort!(intervals_a)
    sort!(intervals_b)
    intersections = Any[]
    i = 1
    j = 1
    while i <= length(intervals_a) && j <= length(intervals_b)
        ai = intervals_a[i]
        bj = intervals_b[j]
        if isless(ai.seqname, bj.seqname) || (ai.seqname == bj.seqname && ai.last < bj.first)
            i += 1
        elseif isless(bj.seqname, ai.seqname) || (ai.seqname == bj.seqname && bj.last < ai.first)
            j += 1
        else
            k = j
            while k <= length(intervals_b) && intervals_b[k].first <= ai.last
                if isoverlapping(ai, intervals_b[k]) && filter(ai, intervals_b[k])
                    push!(intersections, (ai, intervals_b[k]))
                end
                k += 1
            end
            i += 1
        end
    end
    return intersections
end

function simple_coverage(intervals)
    seqlens = Dict{String, Int}()
    for interval in intervals
        if get(seqlens, interval.seqname, -1) < interval.last
            seqlens[interval.seqname] = interval.last
        end
    end

    covarrays = Dict{String, Vector{Int}}()
    for (seqname, seqlen) in seqlens
        covarrays[seqname] = zeros(Int, seqlen)
    end

    for interval in intervals
        arr = covarrays[interval.seqname]
        for i in interval.first:interval.last
            arr[i] += 1
        end
    end

    covintervals = GenomicInterval{UInt32}[]
    for (seqname, arr) in covarrays
        i = j = 1
        while i <= length(arr)
            if arr[i] > 0
                j = i + 1
                while j <= length(arr) && arr[j] == arr[i]
                    j += 1
                end
                push!(covintervals, GenomicInterval{UInt32}(seqname, i, j - 1, STRAND_BOTH, arr[i]))
                i = j
            else
                i += 1
            end
        end
    end

    return covintervals
end



@testset "Strand" begin
    @testset "Constructor" begin
        @test Strand('?') === STRAND_NA
        @test Strand('+') === STRAND_POS
        @test Strand('-') === STRAND_NEG
        @test Strand('.') === STRAND_BOTH
        @test_throws Exception Strand('x')
    end

    @testset "Conversion" begin
        @test convert(Strand, '?') === STRAND_NA
        @test convert(Strand, '+') === STRAND_POS
        @test convert(Strand, '-') === STRAND_NEG
        @test convert(Strand, '.') === STRAND_BOTH

        @test convert(Char, STRAND_NA) === '?'
        @test convert(Char, STRAND_POS) === '+'
        @test convert(Char, STRAND_NEG) === '-'
        @test convert(Char, STRAND_BOTH) === '.'
    end

    @testset "Order" begin
        @test STRAND_NA < STRAND_POS < STRAND_NEG < STRAND_BOTH
    end

    @testset "IO" begin
        # show
        buf = IOBuffer()
        for s in [STRAND_NA, STRAND_POS, STRAND_NEG, STRAND_BOTH]
            show(buf, s); print(buf, " ")
        end
        @test String(take!(buf)) == "STRAND_NA STRAND_POS STRAND_NEG STRAND_BOTH "

        # print
        buf = IOBuffer()
        for s in [STRAND_NA, STRAND_POS, STRAND_NEG, STRAND_BOTH]
            print(buf, s)
        end
        @test String(take!(buf)) == "?+-."

        # read and write
        buf = IOBuffer()
        @test write(buf, STRAND_NA)   == 1
        @test write(buf, STRAND_POS)  == 1
        @test write(buf, STRAND_NEG)  == 1
        @test write(buf, STRAND_BOTH) == 1
        seekstart(buf)
        @test read(buf, Strand) === STRAND_NA
        @test read(buf, Strand) === STRAND_POS
        @test read(buf, Strand) === STRAND_NEG
        @test read(buf, Strand) === STRAND_BOTH
        @test eof(buf)
    end
end

@testset "GenomicInterval" begin
    @testset "Constructor" begin
        i = GenomicInterval("chr1", 10, 20)
        @test seqname(i) == "chr1"
        @test leftposition(i) == 10
        @test rightposition(i) == 20
        @test strand(i) == STRAND_BOTH
        @test i == GenomicInterval("chr1", 10:20)

        i1 = GenomicInterval("chr1", 10, 20, '+')
        i2 = GenomicInterval("chr1", 10, 20, STRAND_POS)
        @test i1 == i2
        @test i1 == GenomicInterval("chr1", 10:20, '+')

        i1 = GenomicInterval("chr2", 5692667, 5701385, '+',        "SOX11")
        i2 = GenomicInterval("chr2", 5692667, 5701385, STRAND_POS, "SOX11")
        @test i1 == i2
        @test i1 == GenomicInterval("chr2", 5692667:5701385, '+', "SOX11")
    end
end

@testset "GenomicIntervalCollection" begin

    @testset "Constructor" begin
        @test GenomicIntervalCollection{Int}() == GenomicIntervalCollection{GenomicInterval{Int}}()

        intervals = [GenomicInterval("test", 1, 2)]
        @test GenomicIntervalCollection{Nothing}(intervals) == GenomicIntervalCollection{GenomicInterval{Nothing}}(intervals)
    end

    @testset "Insertion/Iteration" begin
        n = 100000
        intervals = random_intervals(["one", "two", "three"], 1000000, n)
        ic = GenomicIntervalCollection{Int}()

        @test isempty(ic)
        @test collect(GenomicInterval{Int}, ic) == GenomicInterval{Int}[]

        for interval in intervals
            push!(ic, interval)
        end
        @test is_all_ordered(collect(GenomicInterval{Int}, ic))
    end

    @testset "Intersection" begin
        n = 1000
        Random.seed!(1234)
        intervals_a = random_intervals(["one", "two", "three"], 1000000, n)
        intervals_b = random_intervals(["one", "three", "four"], 1000000, n)

        # empty versus empty
        ic_a = GenomicIntervalCollection{Int}()
        ic_b = GenomicIntervalCollection{Int}()
        @test collect(eachoverlap(ic_a, ic_b)) == Any[]

        # empty versus non-empty
        for interval in intervals_a
            push!(ic_a, interval)
        end

        @test collect(eachoverlap(ic_a, ic_b)) == Any[]
        @test collect(eachoverlap(ic_b, ic_a)) == Any[]

        # non-empty versus non-empty
        for interval in intervals_b
            push!(ic_b, interval)
        end

        @test sort(collect(eachoverlap(ic_a, ic_b))) == sort(simple_intersection(intervals_a, intervals_b))
        @test sort(collect(eachoverlap(ic_a, ic_b, filter=(a,b) -> isodd(first(a))))) == sort(simple_intersection(intervals_a, intervals_b, filter=(a,b) -> isodd(first(a))))
    end

    @testset "Show" begin
        ic = GenomicIntervalCollection{Int}()
        show(devnull, ic)

        push!(ic, GenomicInterval{Int}("one", 1, 1000, STRAND_POS, 0))
        show(devnull, ic)

        intervals = random_intervals(["one", "two", "three"], 1000000, 100)
        for interval in intervals
            push!(ic, interval)
        end
        show(devnull, ic)

        show(devnull, STRAND_NA)
        show(devnull, STRAND_POS)
        show(devnull, STRAND_NEG)
        show(devnull, STRAND_BOTH)
    end
end

@testset "GenomicIntervalStream" begin
    @testset "Intersection" begin
        n = 1000
        Random.seed!(1234)
        intervals_a = random_intervals(["one", "two", "three"], 1000000, n)
        intervals_b = random_intervals(["one", "three", "four"], 1000000, n)

        ic_a = GenomicIntervalCollection{Int}()
        ic_b = GenomicIntervalCollection{Int}()

        for interval in intervals_a
            push!(ic_a, interval)
        end

        for interval in intervals_b
            push!(ic_b, interval)
        end

        # non-empty versus non-empty, stream intersection
        it = eachoverlap(ic_a, ic_b, isless)
        @test sort(collect(it)) == sort(simple_intersection(intervals_a, intervals_b))

        it = eachoverlap(ic_a, ic_b, isless, filter=(a,b) -> isodd(first(a)))
        @test sort(collect(it)) == sort(simple_intersection(intervals_a, intervals_b, filter=(a,b) -> isodd(first(a))))

        it = eachoverlap(
            [GenomicInterval("a", 1, 100, STRAND_POS, nothing), GenomicInterval("c", 1, 100, STRAND_POS, nothing)],
            [GenomicInterval("a", 1, 100, STRAND_POS, nothing), GenomicInterval("b", 1, 100, STRAND_POS, nothing)],
            isless)
        @test length(collect(it)) == 1

        it = eachoverlap(
            [GenomicInterval("c", 1, 100, STRAND_POS, nothing), GenomicInterval("d", 1, 100, STRAND_POS, nothing)],
            [GenomicInterval("b", 1, 100, STRAND_POS, nothing), GenomicInterval("d", 1, 100, STRAND_POS, nothing)],
            isless)
        @test length(collect(it)) == 1

        # unsorted streams are not allowed
        @test_throws Exception begin
            it = eachoverlap(
                [GenomicInterval("b", 1, 1000, STRAND_POS, nothing),
                 GenomicInterval("a", 1, 1000, STRAND_POS, nothing)],
                [GenomicInterval("a", 1, 1000, STRAND_POS, nothing),
                 GenomicInterval("b", 1, 1000, STRAND_POS, nothing)], isless)
            collect(it)
        end

        @test_throws Exception begin
            it = eachoverlap(
                [GenomicInterval("a", 1, 1000, STRAND_POS, nothing),
                 GenomicInterval("a", 500, 1000, STRAND_POS, nothing),
                 GenomicInterval("a", 400, 2000, STRAND_POS, nothing)],
                [GenomicInterval("a", 1, 1000, STRAND_POS, nothing),
                 GenomicInterval("b", 1, 1000, STRAND_POS, nothing)], isless)
            collect(it)
        end
    end

    @testset "GenomicIntervalStream Intersection" begin
        n = 1000
        Random.seed!(1234)
        intervals_a = random_intervals(["one", "two", "three"], 1000000, n)
        intervals_b = random_intervals(["one", "two", "three"], 1000000, n)

        ic_a = GenomicIntervalCollection{Int}()
        ic_b = GenomicIntervalCollection{Int}()

        for interval in intervals_a
            push!(ic_a, interval)
        end

        for interval in intervals_b
            push!(ic_b, interval)
        end

        @test sort(collect(eachoverlap(ic_a, ic_b, isless))) == sort(simple_intersection(intervals_a, intervals_b))
        @test sort(collect(eachoverlap(ic_a, ic_b, isless, filter=(a,b) -> isodd(first(a))))) == sort(simple_intersection(intervals_a, intervals_b, filter=(a,b) -> isodd(first(a))))
    end

    @testset "GenomicIntervalStream Coverage" begin
        n = 10000
        Random.seed!(1234)
        intervals = random_intervals(["one", "two", "three"], 1000000, n)

        ic = GenomicIntervalCollection{Int}()
        for interval in intervals
            push!(ic, interval)
        end

        @test sort(simple_coverage(intervals)) == sort(collect(coverage(ic)))
    end

    @testset "eachoverlap" begin
        i = GenomicInterval("chr1", 1, 10)
        intervals_a = typeof(i)[]
        @test length(collect(eachoverlap(intervals_a, intervals_a))) == 0

        intervals_a = [GenomicInterval("chr1", 1, 10)]
        intervals_b = eltype(intervals_a)[]
        @test length(collect(eachoverlap(intervals_a, intervals_b))) == 0
        @test length(collect(eachoverlap(intervals_b, intervals_a))) == 0

        intervals_a = [GenomicInterval("chr1", 1, 10)]
        @test length(collect(eachoverlap(intervals_a, intervals_a))) == 1

        intervals_a = [GenomicInterval("chr1", 1, 10)]
        intervals_b = [GenomicInterval("chr2", 1, 10)]
        @test length(collect(eachoverlap(intervals_a, intervals_b))) == 0
        @test length(collect(eachoverlap(intervals_b, intervals_a))) == 0

        intervals_a = [GenomicInterval("chr1", 1, 10)]
        intervals_b = [GenomicInterval("chr1", 11, 15)]
        @test length(collect(eachoverlap(intervals_a, intervals_b))) == 0
        @test length(collect(eachoverlap(intervals_b, intervals_a))) == 0

        intervals_a = [GenomicInterval("chr1", 11, 15)]
        intervals_b = [GenomicInterval("chr1", 1, 10), GenomicInterval("chr1", 12, 13)]
        @test length(collect(eachoverlap(intervals_a, intervals_b))) == 1
        @test length(collect(eachoverlap(intervals_b, intervals_a))) == 1

        intervals_a = [GenomicInterval("chr1", 1, 2), GenomicInterval("chr1", 2, 5), GenomicInterval("chr2", 1, 10)]
        intervals_b = [GenomicInterval("chr1", 1, 2), GenomicInterval("chr1", 2, 3), GenomicInterval("chr2", 1, 2)]
        @test length(collect(eachoverlap(intervals_a, intervals_b))) == 5
        @test length(collect(eachoverlap(intervals_b, intervals_a))) == 5

        intervals_a = [GenomicInterval("chr1", 1, 2), GenomicInterval("chr1", 3, 5), GenomicInterval("chr2", 1, 10)]
        @test length(collect(eachoverlap(intervals_a, intervals_a))) == 3

        # compare generic and specific eachoverlap methods
        intervals_a = [GenomicInterval("chr1", 1, 2), GenomicInterval("chr1", 1, 3), GenomicInterval("chr1", 5, 9),
                       GenomicInterval("chr2", 1, 5), GenomicInterval("chr2", 6, 6), GenomicInterval("chr2", 6, 8)]
        intervals_b = intervals_a
        ic_a = GenomicIntervalCollection(intervals_a)
        ic_b = GenomicIntervalCollection(intervals_b)
        iter1 = eachoverlap(intervals_a, intervals_b)
        iter2 = eachoverlap(intervals_a, ic_b)
        iter3 = eachoverlap(ic_a, intervals_b)
        iter4 = eachoverlap(ic_a, ic_b)
        @test collect(iter1) == collect(iter2) == collect(iter3) == collect(iter4)
    end
end

@testset "Custom Concrete Types" begin

    # Define custom Interval type.
    struct GATC{T} <: GenomicFeatures.AbstractGenomicInterval{T}
        seqname::String
        first::Int64
        last::Int64
        metadata::T
    end

    gatcs = [GATC("test1",left,right,nothing) for (left, right) in zip(1:4:3*4, 4:4:3*4)]

    # Collection.
    col_gatc = GenomicIntervalCollection(gatcs[1:2])

    push!(col_gatc, gatcs[3])

    @test collect(col_gatc) == gatcs

    # TODO: Mixed types.
    # push!(gatc_col, GenomicInterval("test1", 9, 12))

    # Overlap.
    intervals_b = intervals_a = gatcs
    ic_a = GenomicIntervalCollection(intervals_a)
    ic_b = GenomicIntervalCollection(intervals_b)
    iter1 = eachoverlap(intervals_a, intervals_b)
    iter2 = eachoverlap(intervals_a, ic_b)
    iter3 = eachoverlap(ic_a, intervals_b)
    iter4 = eachoverlap(ic_a, ic_b)
    @test collect(iter1) == collect(iter2) == collect(iter3) == collect(iter4)

    # Coverage.
    @test collect(coverage(col_gatc)) == [GenomicInterval{UInt32}("test1",1,12,'.',1)] #TODO: relax Number comparisons.

end

@testset "Check Deprecated" begin
    # Interval
    @test (@test_deprecated Interval("test", 1, 2)) == GenomicInterval("test", 1, 2)
    @test (@test_deprecated Interval("test", 1:2)) == GenomicInterval("test", 1, 2)

    @test_deprecated intervals = [Interval("test", 1, 2)]

    # IntervalCollection
    @test (@test_deprecated IntervalCollection{Nothing}()) == GenomicIntervalCollection{Nothing}() == GenomicIntervalCollection{GenomicInterval{Nothing}}()
    @test (@test_deprecated IntervalCollection(intervals)) == (@test_deprecated IntervalCollection{Nothing}(intervals)) == GenomicIntervalCollection(intervals) == GenomicIntervalCollection{Nothing}(intervals) == GenomicIntervalCollection{GenomicInterval{Nothing}}(intervals)

end
