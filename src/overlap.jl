# Overlap Iterator
# ================

struct OverlapIterator{Sa,Sb,F,G}
    intervals_a::Sa
    intervals_b::Sb
    isless::F
    filter::G
end

function Base.eltype(::Type{OverlapIterator{Sa,Sb,F,G}}) where {Sa,Sb,F,G}
    return Tuple{Interval{metadatatype(Sa)},Interval{metadatatype(Sb)}}
end

function Base.IteratorSize(::Type{OverlapIterator{Sa,Sb,F,G}}) where {Sa,Sb,F,G}
    return Base.SizeUnknown()
end

"""
    eachoverlap(intervals_a, intervals_b, [seqname_isless=Base.isless])

Create an iterator of overlapping intervals between `intervals_a` and `intervals_b`.

This function assumes elements of `intervals_a` and `intervals_b` are sorted by
its sequence name and left position.  If the element type is not a subtype of
`GenomicFeatures.Interval`, elements are converted to `Interval` objects.

The third optional argument is a function that defines the order of sequence
names. The default function is `Base.isless`, which is the lexicographical
order.
"""
function eachoverlap(intervals_a, intervals_b, seqname_isless=Base.isless; filter=true_cmp)
    return OverlapIterator(intervals_a, intervals_b, seqname_isless, filter)
end

struct OverlapIteratorState{Sa,Sb,Ta,Tb}
    next_a::Sa
    next_b::Sb
    queue::Queue{Interval{Tb}}
    queue_index::Int
end

function OverlapIteratorState(
        Ta::Type, Tb::Type, next_a::Sa, next_b::Sb, queue::Queue, queue_index::Int) where {Sa, Sb}
    return OverlapIteratorState{Sa,Sb,Ta,Tb}(next_a, next_b, queue, queue_index)
end

function OverlapIteratorState(Ta::Type, Tb::Type, next_a::Sa, next_b::Sb) where {Sa, Sb}
    queue = Queue{Interval{Tb}}()
    return OverlapIteratorState{Sa,Sb,Ta,Tb}(next_a, next_b, queue, 1)
end


function Base.iterate(iter::OverlapIterator)
    next_a = iterate(iter.intervals_a)
    next_b = iterate(iter.intervals_b)

    Ta = metadatatype(iter.intervals_a)
    Tb = metadatatype(iter.intervals_b)
    state = OverlapIteratorState(Ta, Tb, next_a, next_b)

    return iterate(iter, state)
end


# check i1 and i2 are ordered
function check_ordered(i1, i2, compare_func)
    if !isordered(i1, i2, compare_func)
        error("intervals are not sorted")
    end
    return nothing
end


function Base.iterate(iter::OverlapIterator, state::OverlapIteratorState{Sa,Sb,Ta,Tb}) where {Sa,Sb,Ta,Tb}
    next_a      = state.next_a
    next_b      = state.next_b
    queue       = state.queue
    queue_index = state.queue_index

    if next_a === nothing
        return nothing
    end

    entry_a, state_a = next_a
    interval_a = convert(Interval{Ta}, entry_a)

    while true
        if queue_index > lastindex(state.queue)
            # end of queue: add more to queue, or advance a
            if next_b === nothing
                next_a = iterate(iter.intervals_a, state_a)
                if next_a === nothing
                    return break
                end

                entry_a, state_a = next_a
                next_interval_a = convert(Interval{Ta}, entry_a)
                check_ordered(interval_a, next_interval_a, iter.isless)
                interval_a = next_interval_a
                queue_index = firstindex(state.queue)
            else
                entry_b, state_b = next_b
                interval_b = convert(Interval{Tb}, entry_b)
                if !isempty(queue)
                    check_ordered(queue[end], interval_b, iter.isless)
                end
                push!(queue, interval_b)
                next_b = iterate(iter.intervals_b, state_b)
            end
        else
            entry_a, state_a = next_a
            interval_a = convert(Interval{Ta}, entry_a)
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
                next_interval_a = convert(Interval{Ta}, entry_a)

                check_ordered(interval_a, next_interval_a, iter.isless)
                interval_a = next_interval_a
                queue_index = firstindex(state.queue)
            elseif c == 0
                if iter.filter(interval_a, interval_b)
                    return ((interval_a, interval_b),
                        OverlapIteratorState(Ta, Tb, next_a, next_b, queue, queue_index))
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
function compare_overlap(i1::Interval, i2::Interval, isless::Function)
    if isless(i1.seqname, i2.seqname)::Bool
        return -1
    elseif isless(i2.seqname, i1.seqname)::Bool
        return +1
    else  # i1.seqname == i2.seqname
        if i1.last < i2.first
            return -1
        elseif i1.first > i2.last
            return +1
        else
            return 0
        end
    end
end

# Faster comparison for `Base.isless`.  Note that `Base.isless` must be
# consistent wtih `Base.cmp` to work correctly.
function compare_overlap(i1::Interval, i2::Interval, ::typeof(Base.isless))
    c = cmp(i1.seqname, i2.seqname)
    if c != 0
        return c
    end
    if i1.last < i2.first
        return -1
    elseif i1.first > i2.last
        return +1
    else
        return 0
    end
end
