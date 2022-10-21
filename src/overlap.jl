# Overlap Iterator
# ================

struct OverlapIterator{Sa,Sb,L,F}
    interval_stream_a::Sa
    interval_stream_b::Sb
    isless::L
    filter::F
end

function Base.eltype(::Type{OverlapIterator{Sa,Sb,L,F}}) where {Sa,Sb,L,F}
    return Tuple{intervaltype(eltype(Sa)),intervaltype(eltype(Sb))} #Note: The iterator will attempt to convert to these eltypes.
end

function Base.IteratorSize(::Type{OverlapIterator{Sa,Sb,L,F}}) where {Sa,Sb,L,F}
    return Base.SizeUnknown()
end

"""
    eachoverlap(intervals_a, intervals_b, [seqname_isless=Base.isless])

Create an iterator of overlapping intervals between `intervals_a` and `intervals_b`.

This function assumes elements of `intervals_a` and `intervals_b` are sorted by its sequence name and left position.
If the element type is not a subtype of `GenomicFeatures.GenomicInterval`, elements are converted to `GenomicInterval` objects.

The third optional argument is a function that defines the order of sequence names.
The default function is `Base.isless`, which is the lexicographical order.
"""
function eachoverlap(intervals_a, intervals_b, seqname_isless=Base.isless; filter=true_cmp)
    return OverlapIterator(intervals_a, intervals_b, seqname_isless, filter)
end

struct OverlapIteratorState{Na,Nb,Ia,Ib}
    next_a::Na # Note: of the form (value, state)
    next_b::Nb
    queue::Queue{Ib}
    queue_index::Int
end

function OverlapIteratorState(Ia::Type, Ib::Type, next_a::Na, next_b::Nb, queue::Queue, queue_index::Int) where {Na, Nb}
    return OverlapIteratorState{Na,Nb,Ia,Ib}(next_a, next_b, queue, queue_index)
end

function OverlapIteratorState(Ia::Type, Ib::Type, next_a::Na, next_b::Nb) where {Na, Nb}
    queue = Queue{Ib}()
    return OverlapIteratorState{Na,Nb,Ia,Ib}(next_a, next_b, queue, 1)
end


function Base.iterate(iter::OverlapIterator)
    next_a = iterate(iter.interval_stream_a)
    next_b = iterate(iter.interval_stream_b)

    Ia = intervaltype(eltype(iter.interval_stream_a))
    Ib = intervaltype(eltype(iter.interval_stream_b))

    state = OverlapIteratorState(Ia, Ib, next_a, next_b)

    return iterate(iter, state)
end


# check i1 and i2 are ordered
function check_ordered(i1, i2, compare_func)
    if !isordered(i1, i2, compare_func)
        error("intervals are not sorted")
    end
    return nothing
end


function Base.iterate(iter::OverlapIterator, state::OverlapIteratorState{Na,Nb,Ia,Ib}) where {Na,Nb,Ia,Ib}
    next_a      = state.next_a
    next_b      = state.next_b
    queue       = state.queue
    queue_index = state.queue_index

    if next_a === nothing
        return nothing
    end

    entry_a, state_a = next_a
    interval_a = convert(Ia, entry_a)

    while true
        if queue_index > lastindex(state.queue)
            # end of queue: add more to queue, or advance a
            if next_b === nothing
                next_a = iterate(iter.interval_stream_a, state_a)
                if next_a === nothing
                    return break
                end

                entry_a, state_a = next_a
                next_interval_a = convert(Ia, entry_a)
                check_ordered(interval_a, next_interval_a, iter.isless)
                interval_a = next_interval_a
                queue_index = firstindex(state.queue)
            else
                entry_b, state_b = next_b
                interval_b = convert(Ib, entry_b)
                if !isempty(queue)
                    check_ordered(queue[end], interval_b, iter.isless)
                end
                push!(queue, interval_b)
                next_b = iterate(iter.interval_stream_b, state_b)
            end
        else
            entry_a, state_a = next_a
            interval_a = convert(Ia, entry_a)
            interval_b = queue[queue_index]
            c = compare_overlap(interval_a, interval_b, iter.isless)
            queue_index += 1

            if c < 0
                # No more possible intersections with interval_a, advance
                next_a = iterate(iter.interval_stream_a, state_a)
                if next_a === nothing
                    break
                end
                entry_a, state_a = next_a
                next_interval_a = convert(Ia, entry_a)

                check_ordered(interval_a, next_interval_a, iter.isless)
                interval_a = next_interval_a
                queue_index = firstindex(state.queue)
            elseif c == 0
                if iter.filter(interval_a, interval_b)
                    return ((interval_a, interval_b), OverlapIteratorState(Ia, Ib, next_a, next_b, queue, queue_index))
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
function compare_overlap(i1::GenomicInterval, i2::GenomicInterval, isless::Function)
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
function compare_overlap(i1::GenomicInterval, i2::GenomicInterval, ::typeof(Base.isless))
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
