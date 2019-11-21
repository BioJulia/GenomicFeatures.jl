# Overlap Iterator
# ================

struct OverlapIterator{A,B,F,G}
    intervals_a::A # Stream A.
    intervals_b::B # Stream B.
    isless::F
    filter::G
end

function Base.eltype(::Type{OverlapIterator{A,B,F,G}}) where {A,B,F,G}
    return Tuple{GenomicInterval{metadatatype(A)},GenomicInterval{metadatatype(B)}}
end

function Base.eltype(::Type{OverlapIterator{A,B,F,G}}) where {Ia<:AbstractGenomicInterval, A<:Union{AbstractVector{Ia},GenomicIntervalCollection{Ia}}, Ib<:AbstractGenomicInterval,B<:Union{AbstractVector{Ib},GenomicIntervalCollection{Ib}},F,G}
    return Tuple{Ia,Ib}
end

function Base.IteratorSize(::Type{OverlapIterator{A,B,F,G}}) where {A,B,F,G}
    return Base.SizeUnknown()
end

"""
    eachoverlap(intervals_a, intervals_b, [seqname_isless=Base.isless])

Create an iterator of overlapping intervals between `intervals_a` and `intervals_b`.

This function assumes elements of `intervals_a` and `intervals_b` are sorted by its sequence name and left position.
If the element type is not a subtype of `GenomicFeatures.AbstractGenomicInterval`, elements are converted to `GenomicInterval` objects.

The third optional argument is a function that defines the order of sequence names.
The default function is `Base.isless`, which is the lexicographical order.
"""
function eachoverlap(intervals_a, intervals_b, seqname_isless=Base.isless; filter=true_cmp)
    return OverlapIterator(intervals_a, intervals_b, seqname_isless, filter)
end

struct OverlapIteratorState{Na,Nb,Ea,Eb}
    next_a::Na #Note: Na <: Union{Nothing, Tuple{Ea, Sa}}, where Ea <: AbstractGenomicInterval, and state Sa <: Any
    next_b::Nb #Note: Nb <: Union{Nothing, Tuple{Eb, Sb}}, where Eb <: AbstractGenomicInterval, and state Sb <: Any
    queue::Queue{Union{Ea,Eb}}
    queue_index::Int
end

function OverlapIteratorState{Na,Nb,Ea,Eb}(next_a::Na, next_b::Nb) where {Na,Nb,Ea,Eb}
    Ia = GenomicInterval{Ea}
    Ib = GenomicInterval{Eb}

    queue = Queue{Union{Ia,Ib}}()
    return OverlapIteratorState{Na,Nb,Ea,Eb}(next_a, next_b, queue, 1)
end

function OverlapIteratorState{Na,Nb,Ea,Eb}(next_a::Na, next_b::Nb) where {Na,Nb,Ea<:AbstractGenomicInterval,Eb<:AbstractGenomicInterval}
    queue = Queue{Union{Ea,Eb}}()
    return OverlapIteratorState{Na,Nb,Ea,Eb}(next_a, next_b, queue, 1)
end

function OverlapIteratorState(Ea::Type, Eb::Type, next_a::Na, next_b::Nb) where {Na, Nb}
    return OverlapIteratorState{Na,Nb,Ea,Eb}(next_a, next_b)
end

function OverlapIteratorState(Ea::Type, Eb::Type, next_a::Na, next_b::Nb, queue::Queue, queue_index::Int) where {Na, Nb}
    return OverlapIteratorState{Na,Nb,Ea,Eb}(next_a, next_b, queue, queue_index)
end

# Initial iteration.
function Base.iterate(iter::OverlapIterator{A,B,F,G}) where {A,B,F,G}
    next_a = iterate(iter.intervals_a) #Note: returns (value, state) or nothing.
    next_b = iterate(iter.intervals_b)

    Ea = eltype(iter.intervals_a) #TODO: use eltype of OverlapIterator?
    Eb = eltype(iter.intervals_b)

    state = OverlapIteratorState(Ea, Eb, next_a, next_b) #TODO: consider doing next_a and next-b conversion here.

    return iterate(iter, state)
end


# check i1 and i2 are ordered
function check_ordered(i1, i2, compare_func)
    if !isordered(i1, i2, compare_func)
        error("intervals are not sorted")
    end
    return nothing
end

# Subsequent iteration.
function Base.iterate(iter::OverlapIterator{A,B,F,G}, state::OverlapIteratorState{Na,Nb,Ea,Eb}) where {A,B,F,G,Na,Nb,Ea,Eb}
    next_a      = state.next_a
    next_b      = state.next_b
    queue       = state.queue
    queue_index = state.queue_index

    if next_a === nothing
        return nothing
    end

    entry_a, state_a = next_a
    interval_a = typeof(entry_a) <: Ea ? entry_a : convert(Ea, entry_a) #TODO: handle conversion elsewhere.

    while true
        if queue_index > lastindex(queue)
            # end of queue: add more to queue, or advance a
            if next_b === nothing
                next_a = iterate(iter.intervals_a, state_a)
                if next_a === nothing
                    return break
                end

                entry_a, state_a = next_a
                next_interval_a = typeof(entry_a) <: Ea ? entry_a : convert(Ea, entry_a) #TODO: handle conversion elsewhere.
                check_ordered(interval_a, next_interval_a, iter.isless)
                interval_a = next_interval_a
                queue_index = firstindex(queue)
            else
                entry_b, state_b = next_b
                interval_b = typeof(entry_b) <: Eb ? entry_b : convert(Eb, entry_b) #TODO: handle conversion elsewhere.
                if !isempty(queue)
                    check_ordered(queue[end], interval_b, iter.isless)
                end
                push!(queue, interval_b)
                next_b = iterate(iter.intervals_b, state_b)
            end
        else
            entry_a, state_a = next_a
            interval_a = typeof(entry_a) <: Ea ? entry_a : convert(Ea, entry_a) #TODO: handle conversion elsewhere.
            interval_b = queue[queue_index]
            c = compare_overlap(interval_a, interval_b, iter.isless)
            queue_index += 1

            if c < 0
                # No more possible intersections with interval_a, advance
                next_a = iterate(iter.intervals_a, state_a)
                if next_a === nothing
                    break
                end
                entry_a, state_a = next_a
                next_interval_a = typeof(entry_a) <: Ea ? entry_a : convert(Ea, entry_a) #TODO: handle conversion elsewhere.

                check_ordered(interval_a, next_interval_a, iter.isless)
                interval_a = next_interval_a
                queue_index = firstindex(queue)
            elseif c == 0
                if iter.filter(interval_a, interval_b)
                    return ((interval_a, interval_b), OverlapIteratorState(Ea, Eb, next_a, next_b, queue, queue_index))
                end
            else
                if queue_index == firstindex(queue) + 1
                    # noting else can intersect front of queue
                    popfirst!(queue)
                end
            end
        end
    end

    # no more intersections found
    return nothing
end

# Return:
#   -1 when `i1` precedes `i2`,
#   0 when `i1` overlaps with `i2`, and
#   +1 when `i1` follows `i2`.
function compare_overlap(i1::AbstractGenomicInterval, i2::AbstractGenomicInterval, isless::Function)
    if isless(seqname(i1), seqname(i2))
        return -1
    end

    if isless(seqname(i2), seqname(i1))
        return +1
    end

    # seqname(i1) == seqname(i2)
    if rightposition(i1) < leftposition(i2)
        return -1
    end

    if leftposition(i1) > rightposition(i2)
        return +1
    end

    return 0
end

# Faster comparison for `Base.isless`.  Note that `Base.isless` must be consistent wtih `Base.cmp` to work correctly.
function compare_overlap(i1::AbstractGenomicInterval, i2::AbstractGenomicInterval, ::typeof(Base.isless))
    c = cmp(seqname(i1), seqname(i2))

    if c != 0
        return c
    end

    if rightposition(i1) < leftposition(i2)
        return -1
    end

    if leftposition(i1) > rightposition(i2)
        return +1
    end

    return 0
end
