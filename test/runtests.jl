using GenomicFeatures
using Test

using Distributions
import Random


# Test that an array of intervals is well ordered
function is_all_ordered(intervals::Vector{I}) where I <: Interval
    for i = 2:length(intervals)
        if !GenomicFeatures.isordered(intervals[i-1], intervals[i])
            return false
        end
    end
    return true
end

# Generate an array of n random Interval{Int} object. With sequence names
# samples from seqnames, and intervals drawn to lie in [1, maxpos].
function random_intervals(seqnames, maxpos::Int, n::Int)
    seq_dist = Categorical(length(seqnames))
    strand_dist = Categorical(2)
    length_dist = Normal(1000, 1000)
    intervals = Vector{Interval{Int}}(undef, n)
    for i in 1:n
        intlen = maxpos
        while intlen >= maxpos || intlen <= 0
            intlen = ceil(Int, rand(length_dist))
        end
        first = rand(1:maxpos-intlen)
        last = first + intlen - 1
        strand = rand(strand_dist) == 1 ? STRAND_POS : STRAND_NEG
        intervals[i] = Interval{Int}(seqnames[rand(seq_dist)],
                                     first, last, strand, i)
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
        if isless(ai.seqname, bj.seqname) ||
           (ai.seqname == bj.seqname && ai.last < bj.first)
            i += 1
        elseif isless(bj.seqname, ai.seqname) ||
               (ai.seqname == bj.seqname && bj.last < ai.first)
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

    covintervals = Interval{UInt32}[]
    for (seqname, arr) in covarrays
        i = j = 1
        while i <= length(arr)
            if arr[i] > 0
                j = i + 1
                while j <= length(arr) && arr[j] == arr[i]
                    j += 1
                end
                push!(covintervals,
                      Interval{UInt32}(seqname, i, j - 1, STRAND_BOTH, arr[i]))
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

@testset "Interval" begin
    @testset "Constructor" begin
        i = Interval("chr1", 10, 20)
        @test seqname(i) == "chr1"
        @test leftposition(i) == 10
        @test rightposition(i) == 20
        @test strand(i) == STRAND_BOTH
        @test i == Interval("chr1", 10:20)

        i1 = Interval("chr1", 10, 20, '+')
        i2 = Interval("chr1", 10, 20, STRAND_POS)
        @test i1 == i2
        @test i1 == Interval("chr1", 10:20, '+')

        i1 = Interval("chr2", 5692667, 5701385, '+',        "SOX11")
        i2 = Interval("chr2", 5692667, 5701385, STRAND_POS, "SOX11")
        @test i1 == i2
        @test i1 == Interval("chr2", 5692667:5701385, '+', "SOX11")
    end
end

@testset "IntervalCollection" begin
    @testset "Insertion/Iteration" begin
        n = 100000
        intervals = random_intervals(["one", "two", "three"], 1000000, n)
        ic = IntervalCollection{Int}()

        @test isempty(ic)
        @test collect(Interval{Int}, ic) == Interval{Int}[]

        for interval in intervals
            push!(ic, interval)
        end
        @test is_all_ordered(collect(Interval{Int}, ic))
    end

    @testset "Intersection" begin
        n = 1000
        Random.seed!(1234)
        intervals_a = random_intervals(["one", "two", "three"], 1000000, n)
        intervals_b = random_intervals(["one", "three", "four"], 1000000, n)

        # empty versus empty
        ic_a = IntervalCollection{Int}()
        ic_b = IntervalCollection{Int}()
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

        @test sort(collect(eachoverlap(ic_a, ic_b))) ==
              sort(simple_intersection(intervals_a, intervals_b))
        @test sort(collect(eachoverlap(ic_a, ic_b, filter=(a,b) -> isodd(first(a))))) ==
              sort(simple_intersection(intervals_a, intervals_b, filter=(a,b) -> isodd(first(a))))
    end

    @testset "Show" begin
        ic = IntervalCollection{Int}()
        show(devnull, ic)

        push!(ic, Interval{Int}("one", 1, 1000, STRAND_POS, 0))
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

@testset "IntervalStream" begin
    @testset "Intersection" begin
        n = 1000
        Random.seed!(1234)
        intervals_a = random_intervals(["one", "two", "three"], 1000000, n)
        intervals_b = random_intervals(["one", "three", "four"], 1000000, n)

        ic_a = IntervalCollection{Int}()
        ic_b = IntervalCollection{Int}()

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
        @test sort(collect(it)) ==
              sort(simple_intersection(intervals_a, intervals_b, filter=(a,b) -> isodd(first(a))))

        it = eachoverlap(
            [Interval("a", 1, 100, STRAND_POS, nothing), Interval("c", 1, 100, STRAND_POS, nothing)],
            [Interval("a", 1, 100, STRAND_POS, nothing), Interval("b", 1, 100, STRAND_POS, nothing)],
            isless)
        @test length(collect(it)) == 1

        it = eachoverlap(
            [Interval("c", 1, 100, STRAND_POS, nothing), Interval("d", 1, 100, STRAND_POS, nothing)],
            [Interval("b", 1, 100, STRAND_POS, nothing), Interval("d", 1, 100, STRAND_POS, nothing)],
            isless)
        @test length(collect(it)) == 1

        # unsorted streams are not allowed
        @test_throws Exception begin
            it = eachoverlap(
                [Interval("b", 1, 1000, STRAND_POS, nothing),
                 Interval("a", 1, 1000, STRAND_POS, nothing)],
                [Interval("a", 1, 1000, STRAND_POS, nothing),
                 Interval("b", 1, 1000, STRAND_POS, nothing)], isless)
            collect(it)
        end

        @test_throws Exception begin
            it = eachoverlap(
                [Interval("a", 1, 1000, STRAND_POS, nothing),
                 Interval("a", 500, 1000, STRAND_POS, nothing),
                 Interval("a", 400, 2000, STRAND_POS, nothing)],
                [Interval("a", 1, 1000, STRAND_POS, nothing),
                 Interval("b", 1, 1000, STRAND_POS, nothing)], isless)
            collect(it)
        end
    end

    @testset "IntervalStream Intersection" begin
        n = 1000
        Random.seed!(1234)
        intervals_a = random_intervals(["one", "two", "three"], 1000000, n)
        intervals_b = random_intervals(["one", "two", "three"], 1000000, n)

        ic_a = IntervalCollection{Int}()
        ic_b = IntervalCollection{Int}()

        for interval in intervals_a
            push!(ic_a, interval)
        end

        for interval in intervals_b
            push!(ic_b, interval)
        end

        @test sort(collect(eachoverlap(ic_a, ic_b, isless))) ==
              sort(simple_intersection(intervals_a, intervals_b))
        @test sort(collect(eachoverlap(ic_a, ic_b, isless, filter=(a,b) -> isodd(first(a))))) ==
              sort(simple_intersection(intervals_a, intervals_b, filter=(a,b) -> isodd(first(a))))
    end

    @testset "IntervalStream Coverage" begin
        n = 10000
        Random.seed!(1234)
        intervals = random_intervals(["one", "two", "three"], 1000000, n)

        ic = IntervalCollection{Int}()
        for interval in intervals
            push!(ic, interval)
        end

        @test sort(simple_coverage(intervals)) == sort(collect(coverage(ic)))
    end

    @testset "eachoverlap" begin
        i = Interval("chr1", 1, 10)
        intervals_a = typeof(i)[]
        @test length(collect(eachoverlap(intervals_a, intervals_a))) == 0

        intervals_a = [Interval("chr1", 1, 10)]
        intervals_b = eltype(intervals_a)[]
        @test length(collect(eachoverlap(intervals_a, intervals_b))) == 0
        @test length(collect(eachoverlap(intervals_b, intervals_a))) == 0

        intervals_a = [Interval("chr1", 1, 10)]
        @test length(collect(eachoverlap(intervals_a, intervals_a))) == 1

        intervals_a = [Interval("chr1", 1, 10)]
        intervals_b = [Interval("chr2", 1, 10)]
        @test length(collect(eachoverlap(intervals_a, intervals_b))) == 0
        @test length(collect(eachoverlap(intervals_b, intervals_a))) == 0

        intervals_a = [Interval("chr1", 1, 10)]
        intervals_b = [Interval("chr1", 11, 15)]
        @test length(collect(eachoverlap(intervals_a, intervals_b))) == 0
        @test length(collect(eachoverlap(intervals_b, intervals_a))) == 0

        intervals_a = [Interval("chr1", 11, 15)]
        intervals_b = [Interval("chr1", 1, 10), Interval("chr1", 12, 13)]
        @test length(collect(eachoverlap(intervals_a, intervals_b))) == 1
        @test length(collect(eachoverlap(intervals_b, intervals_a))) == 1

        intervals_a = [Interval("chr1", 1, 2), Interval("chr1", 2, 5), Interval("chr2", 1, 10)]
        intervals_b = [Interval("chr1", 1, 2), Interval("chr1", 2, 3), Interval("chr2", 1, 2)]
        @test length(collect(eachoverlap(intervals_a, intervals_b))) == 5
        @test length(collect(eachoverlap(intervals_b, intervals_a))) == 5

        intervals_a = [Interval("chr1", 1, 2), Interval("chr1", 3, 5), Interval("chr2", 1, 10)]
        @test length(collect(eachoverlap(intervals_a, intervals_a))) == 3

        # compare generic and specific eachoverlap methods
        intervals_a = [Interval("chr1", 1, 2), Interval("chr1", 1, 3), Interval("chr1", 5, 9),
                       Interval("chr2", 1, 5), Interval("chr2", 6, 6), Interval("chr2", 6, 8)]
        intervals_b = intervals_a
        ic_a = IntervalCollection(intervals_a)
        ic_b = IntervalCollection(intervals_b)
        iter1 = eachoverlap(intervals_a, intervals_b)
        iter2 = eachoverlap(intervals_a, ic_b)
        iter3 = eachoverlap(ic_a, intervals_b)
        iter4 = eachoverlap(ic_a, ic_b)
        @test collect(iter1) == collect(iter2) == collect(iter3) == collect(iter4)
    end
end
