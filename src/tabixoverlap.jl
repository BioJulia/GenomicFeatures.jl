# Tabix Overlap Iterator
# ======================

immutable TabixOverlapIterator{T}
    reader::T
    interval::Interval
end

function Base.eltype{T}(::Type{TabixOverlapIterator{T}})
    return eltype(T)
end

function Base.iteratorsize{T}(::Type{TabixOverlapIterator{T}})
    return Base.SizeUnknown()
end

function GenomicFeatures.eachoverlap(reader::Union{BED.Reader,GFF3.Reader}, interval::Interval)
    if isnull(reader.index)
        throw(ArgumentError("index is null"))
    end
    return TabixOverlapIterator(reader, interval)
end

type TabixOverlapIteratorState{T}
    chunks::Vector{Indexes.Chunk}
    chunkid::Int
    done::Bool
    record::T
end

function Base.start(iter::TabixOverlapIterator)
    @assert !isnull(iter.reader.index)
    return TabixOverlapIteratorState(
        Indexes.overlapchunks(get(iter.reader.index), iter.interval),
        0,
        false,
        eltype(iter)())
end

function Base.done(iter::TabixOverlapIterator, state)
    source = BioCore.IO.stream(iter.reader).source
    if state.chunkid == 0 && !isempty(state.chunks)
        state.chunkid += 1
        seek(source, state.chunks[state.chunkid].start)
    end
    while state.chunkid ≤ endof(state.chunks)
        chunk = state.chunks[state.chunkid]
        while BGZFStreams.virtualoffset(source) < chunk.stop
            read!(iter.reader, state.record)
            c = icmp(state.record, iter.interval)
            if c == 0  # overlapping
                return false
            elseif c > 0
                # no more overlapping records in this chunk
                break
            end
        end
        state.chunkid += 1
        if state.chunkid ≤ endof(state.chunks)
            seek(source, state.chunks[state.chunkid].start)
        end
    end
    # no more overlapping records
    return true
end

function Base.next(iter::TabixOverlapIterator, state)
    return copy(state.record), state
end

function icmp(record, interval)
    c = cmp(seqname(record), interval.seqname)
    if c < 0 || (c == 0 && rightposition(record) < interval.first)
        return -1
    elseif c > 0 || (c == 0 && leftposition(record) > interval.last)
        return +1
    else
        return 0
    end
end
