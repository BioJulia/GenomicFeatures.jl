# Overlap Iterator
# ================

struct OverlapIterator{Sa,Sb,F,G}
    intervals_a::Sa
    intervals_b::Sb
    isless::F
    filter::G
end

function Base.eltype(::Type{OverlapIterator{Sa,Sb,F,G}}) where {Sa,Sb,F,G}
    return Tuple{GenomicInterval{metadatatype(Sa)},GenomicInterval{metadatatype(Sb)}}
end

function Base.eltype(::Type{OverlapIterator{Sa,Sb,F,G}}) where {Ia<:AbstractGenomicInterval, Sa<:Union{AbstractVector{Ia},GenomicIntervalCollection{Ia}}, Ib<:AbstractGenomicInterval,Sb<:Union{AbstractVector{Ib},GenomicIntervalCollection{Ib}},F,G}
    return Tuple{Ia,Ib}
end

function Base.IteratorSize(::Type{OverlapIterator{Sa,Sb,F,G}}) where {Sa,Sb,F,G}
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

struct OverlapIteratorState{Sa,Sb,Ta,Tb}
    next_a::Sa
    next_b::Sb
    queue::Queue{<:AbstractGenomicInterval{Tb}}
    queue_index::Int
end

function OverlapIteratorState(Ta::Type, Tb::Type, next_a::Sa, next_b::Sb, queue::Queue, queue_index::Int) where {Sa, Sb}
    return OverlapIteratorState{Sa,Sb,Ta,Tb}(next_a, next_b, queue, queue_index)
end

function OverlapIteratorState(Ta::Type, Tb::Type, next_a::Sa, next_b::Sb) where {Sa, Sb}
    queue = Queue{GenomicInterval{Tb}}()
    return OverlapIteratorState{Sa,Sb,Ta,Tb}(next_a, next_b, queue, 1)
end

function OverlapIteratorState(Ta::Type, Tb::Type, next_a::Sa, next_b::Sb) where {Sa, Ib<:AbstractGenomicInterval, Sb<:Tuple{Ib, Number}}
    queue = Queue{Ib}()
    return OverlapIteratorState{Sa,Sb,Ta,Tb}(next_a, next_b, queue, 1)
end

function Base.iterate(iter::OverlapIterator)
    next_a = iterate(iter.intervals_a)
    next_b = iterate(iter.intervals_b)

    Ta = metadatatype(iter.intervals_a)
    Tb = metadatatype(iter.intervals_b)

    state = OverlapIteratorState(Ta, Tb, next_a, next_b) #TODO: consider doing next_a and next-b conversion here.

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
    interval_a = typeof(entry_a) <: AbstractGenomicInterval{Ta} ? entry_a : convert(GenomicInterval{Ta}, entry_a) #TODO: use AbstractGenomicInterval and retrieve interval type.

    while true
        if queue_index > lastindex(state.queue)
            # end of queue: add more to queue, or advance a
            if next_b === nothing
                next_a = iterate(iter.intervals_a, state_a)
                if next_a === nothing
                    return break
                end

                entry_a, state_a = next_a
                next_interval_a = typeof(entry_a) <: AbstractGenomicInterval{Ta} ? entry_a : convert(GenomicInterval{Ta}, entry_a) #TODO: use AbstractGenomicInterval and retrieve interval type.
                check_ordered(interval_a, next_interval_a, iter.isless)
                interval_a = next_interval_a
                queue_index = firstindex(state.queue)
            else
                entry_b, state_b = next_b
                interval_b = typeof(entry_b) <: AbstractGenomicInterval{Tb} ? entry_b : convert(GenomicInterval{Tb}, entry_b) #TODO: use AbstractGenomicInterval and retrieve interval type.
                if !isempty(queue)
                    check_ordered(queue[end], interval_b, iter.isless)
                end
                push!(queue, interval_b)
                next_b = iterate(iter.intervals_b, state_b)
            end
        else
            entry_a, state_a = next_a
            interval_a = typeof(entry_a) <: AbstractGenomicInterval{Ta} ? entry_a : convert(GenomicInterval{Ta}, entry_a) #TODO: use AbstractGenomicInterval and retrieve interval type.
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
                next_interval_a = typeof(entry_a) <: AbstractGenomicInterval{Ta} ? entry_a : convert(GenomicInterval{Ta}, entry_a) #TODO: use AbstractGenomicInterval and retrieve interval type.

                check_ordered(interval_a, next_interval_a, iter.isless)
                interval_a = next_interval_a
                queue_index = firstindex(state.queue)
            elseif c == 0
                if iter.filter(interval_a, interval_b)
                    return ((interval_a, interval_b), OverlapIteratorState(Ta, Tb, next_a, next_b, queue, queue_index))
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
    if isless(seqname(i1), seqname(i2))::Bool
        return -1
    elseif isless(seqname(i2), seqname(i1))::Bool
        return +1
    else  # seqname(i1) == seqname(i2)
        if rightposition(i1) < leftposition(i2)
            return -1
        elseif leftposition(i1) > rightposition(i2)
            return +1
        else
            return 0
        end
    end
end

# Faster comparison for `Base.isless`.  Note that `Base.isless` must be consistent wtih `Base.cmp` to work correctly.
function compare_overlap(i1::AbstractGenomicInterval, i2::AbstractGenomicInterval, ::typeof(Base.isless))
    c = cmp(seqname(i1), seqname(i2))
    if c != 0
        return c
    end
    if rightposition(i1) < leftposition(i2)
        return -1
    elseif leftposition(i1) > rightposition(i2)
        return +1
    else
        return 0
    end
end
