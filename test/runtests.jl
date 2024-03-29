using GenomicFeatures
using Test

import Random

include(joinpath(@__DIR__, "Utilities.jl"))
import .Utilities: random_intervals


# Test that an array of intervals is well ordered
function is_all_ordered(intervals::Vector{I}) where I <: GenomicInterval
    for i = 2:length(intervals)
        if !GenomicFeatures.isordered(intervals[i-1], intervals[i])
            return false
        end
    end
    return true
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
        if isless(groupname(ai), groupname(bj)) || (groupname(ai) == groupname(bj) && rightposition(ai) < leftposition(bj))
            i += 1
        elseif isless(groupname(bj), groupname(ai)) || (groupname(ai) == groupname(bj) && rightposition(bj) < leftposition(ai))
            j += 1
        else
            k = j
            while k <= length(intervals_b) && leftposition(intervals_b[k]) <= rightposition(ai)
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
        if get(seqlens, groupname(interval), -1) < rightposition(interval)
            seqlens[groupname(interval)] = rightposition(interval)
        end
    end

    covarrays = Dict{String, Vector{Int}}()
    for (groupname, seqlen) in seqlens
        covarrays[groupname] = zeros(Int, seqlen)
    end

    for interval in intervals
        arr = covarrays[groupname(interval)]
        for i in leftposition(interval):rightposition(interval)
            arr[i] += 1
        end
    end

    covintervals = GenomicInterval{UInt32}[]
    for (groupname, arr) in covarrays
        i = j = 1
        while i <= length(arr)
            if arr[i] > 0
                j = i + 1
                while j <= length(arr) && arr[j] == arr[i]
                    j += 1
                end
                push!(covintervals, GenomicInterval{UInt32}(groupname, i, j - 1, STRAND_BOTH, arr[i]))
                i = j
            else
                i += 1
            end
        end
    end

    return covintervals
end

# Used to check show of AbstractGenomicInterval{Nothing}.
struct WithoutMetadatata <: GenomicFeatures.AbstractGenomicInterval{Nothing}
    groupname::String
    first::Int64
    last::Int64
end

@testset "GenomicFeatures" begin

@testset "Strand" begin
    @testset "Constructor" begin
        @test Strand(0b000) === STRAND_NA
        @test Strand(0b001) === STRAND_POS
        @test Strand(0b010) === STRAND_NEG
        @test Strand(0b011) === STRAND_BOTH

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
        @test groupname(i) == "chr1"
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

    @testset "intervaltype" begin

        intervals = GenomicInterval.("test", [1, 2], [1, 2])
        @test intervals |> eltype |> GenomicFeatures.intervaltype == GenomicInterval{Nothing}

        intervals = GenomicInterval.("test", [1, 2], [1, 2], '.', [1.0, 2.0])
        @test intervals |> eltype |> GenomicFeatures.intervaltype == GenomicInterval{Float64}

        # Example `intervaltype` override.
        function GenomicFeatures.intervaltype(::Type{NamedTuple{(:chrom, :left, :right, :value), Tuple{String, Int, Int, T}}}) where T
            return GenomicInterval{T}
        end

        intervals = NamedTuple.(zip(
            :chrom .=> Iterators.repeated("test", 2),
            :left .=> [1, 2],
            :right .=> [1, 2],
            :value .=> [1.0, 2.0],
        ))

        @test intervals |> eltype |> GenomicFeatures.intervaltype == GenomicInterval{Float64}

    end

    @testset "baseintervaltype" begin
        interval = GenomicInterval("test", 1, 2)
        @test typeof(interval) == GenomicInterval{Nothing}
        @test GenomicFeatures.baseintervaltype(interval) === GenomicInterval

        interval = GenomicInterval("test", 1, 2, '.', 3.0)
        @test typeof(interval) == GenomicInterval{Float64}
        @test GenomicFeatures.baseintervaltype(interval) === GenomicInterval
    end

    @testset "metadatatype" begin
        interval = GenomicInterval("test", 1, 2)
        @test typeof(interval) == GenomicInterval{Nothing}
        @test GenomicFeatures.metadatatype(interval) === Nothing

        interval = GenomicInterval("test", 1, 2, '.', 3.0)
        @test typeof(interval) == GenomicInterval{Float64}
        @test GenomicFeatures.metadatatype(interval) === Float64
    end

    @testset "precedes" begin
        interval_a = GenomicInterval("test", 1, 2)
        interval_b = GenomicInterval("test", 2, 3)
        interval_c = GenomicInterval("test", 3, 4)
        interval_d = GenomicInterval("test2", 3, 4)

        @test GenomicFeatures.precedes(interval_a, interval_a) == false
        @test GenomicFeatures.precedes(interval_a, interval_b) == false
        @test GenomicFeatures.precedes(interval_a, interval_c) == true
        @test GenomicFeatures.precedes(interval_c, interval_a) == false
        @test GenomicFeatures.precedes(interval_a, interval_d) == true
    end

    @test span(GenomicInterval("test", 1, 9)) == length(1:9)
    @test GenomicFeatures.volume(GenomicInterval("test", 1, 9, '?', 2.5)) == length(1:9) * 2.5

end

@testset "GenomicPosition" begin
    @testset "Constructor" begin

        p = GenomicPosition("chr1", 1)
        @test groupname(p) == "chr1"
        @test leftposition(p) == 1
        @test rightposition(p) == 1
        @test position(p) == 1

        @test GenomicPosition("chr1", 1) == GenomicPosition("chr1", 1:1)

        @test typeof(GenomicFeatures.metadata(GenomicPosition{Float64}("chr1", 1, 9))) == Float64
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

    @testset "Bulk Insertion" begin
        i = GenomicInterval("chr1", 10, 20, '+', 1)

        # Test equivalency of shorthand and longhand types.
        @test GenomicIntervalCollection([i]) == GenomicIntervalCollection{Int}([i]) == GenomicIntervalCollection{GenomicInterval{Int}}([i])

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

        # Check hasintersection method.
        @test hasintersection(GenomicInterval("test", 1, 1), GenomicIntervalCollection([GenomicInterval("test", 1,2)])) == true
        @test hasintersection(GenomicInterval("test", 1, 1), GenomicIntervalCollection([GenomicInterval("test", 2,2)])) == false

        # Check hasintersection currying.
        @test GenomicInterval("test", 1, 1) |> hasintersection(GenomicIntervalCollection([GenomicInterval("test", 1,2)])) == true
        @test GenomicInterval("test", 1, 1) |> hasintersection(GenomicIntervalCollection([GenomicInterval("test", 2,2)])) == false
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

@testset "Constructor Conversions" begin

    i = GenomicInterval("chr1", 1, 2, '?', 3)
    nt = (chrom="chr1", left=1, right=2, value=3)

    function Base.convert(::Type{GenomicInterval{Int}}, nt::NamedTuple)
        return GenomicInterval{Int}(nt.chrom, nt[2], nt[3], '?', nt[4])
    end

    @test i == GenomicInterval{Int}(nt)

    @test GenomicIntervalCollection([i]) == GenomicIntervalCollection{Int}([nt])

end #testset Constructor Conversions

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

        # Check coverage operation with different Integer type.
        cov = GenomicIntervalCollection{GenomicInterval{Int64}}()
        @test sort(simple_coverage(intervals)) == collect(GenomicFeatures.coverage!(cov, ic))
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

        # non-overlapping query
        @test length(collect(eachoverlap(ic_a, GenomicInterval("X", 0, 0)))) == 0
    end
end

@testset "Mixed types" begin

    i = GenomicInterval("chr1", 1,4)
    p1 = GenomicPosition("chr1", 2)
    p2 = GenomicPosition("chr1", 5)
    p3 = GenomicPosition("chr2", 5)

    # Sorting.
    @test [i, p1, p2, p3] == sort([p2, p1, p3, i])

    # Push out of order mixed types.
    col = GenomicIntervalCollection{GenomicFeatures.AbstractGenomicInterval{Nothing}}()
    push!(col, p2)
    push!(col, p1)
    push!(col, p3)
    push!(col, i)

    @test collect(col) == [i, p1, p2, p3]

    # Bulk insertion of mixed types.
    col = GenomicIntervalCollection([i, p1, p2, p3])

    @test 4 == length(col)

    # Check overlap.
    interval_p1_equiv = GenomicInterval("chr1", 2,2)
    @test true == isoverlapping(i,interval_p1_equiv)
    @test interval_p1_equiv == p1

    @test true == isoverlapping(i, p1)
    @test false == isoverlapping(i, p2) == isoverlapping(i, p3)

    # Check eachoverlap.
    # @test [(i, p1)] == collect(eachoverlap(i, [p1,p2,p3])) #Note: attempts to broadcast i.
    @test [(i, p1)] == collect(eachoverlap([i], [p1,p2,p3]))
    @test [(i, p1)] == collect(eachoverlap(GenomicIntervalCollection([i]), GenomicIntervalCollection([p1,p2,p3])))
    # @test [(p1, i)] == collect(eachoverlap([p1,p2,p3], i)) #Note: attempts to broadcast i.
    @test [(p1, i)] == collect(eachoverlap([p1,p2,p3], [i]))
    @test [(p1, i)] == collect(eachoverlap(GenomicIntervalCollection([p1,p2,p3]), GenomicIntervalCollection([i])))

    # Pushing mixed types to Queue.
    @test [(i, p1), (p1,p1)] == collect(eachoverlap([i, p1], [p1,p2,p3]))
    @test [(i, p1), (p1,p1)] == collect(eachoverlap(GenomicIntervalCollection([i,p1]), GenomicIntervalCollection([p1,p2,p3])))
    @test [(p1, i), (p1,p1)] == collect(eachoverlap([p1,p2,p3], [i, p1]))
    @test [(p1, i), (p1,p1)] == collect(eachoverlap(GenomicIntervalCollection([p1,p2,p3]), GenomicIntervalCollection([i,p1])))

    # Check coverage.
    @test [GenomicInterval{UInt32}("chr1",1,1,'.',1), GenomicInterval{UInt32}("chr1",2,2,'.',2), GenomicInterval{UInt32}("chr1",3,4,'.',1)] == collect(coverage(GenomicIntervalCollection([i,p1]))) #TODO: relax comparisons.
end


@testset "Custom Concrete Types" begin



    buf = IOBuffer()

    show(buf, WithoutMetadatata("chr1", 1, 2))
    @test String(take!(buf)) == "WithoutMetadatata:\n  group name: chr1\n  leftmost position: 1\n  rightmost position: 2\n  metadata: nothing"

end #testset "Custom Concrete Types"

end #testset GenomicFeatures
