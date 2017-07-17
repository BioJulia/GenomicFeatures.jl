# GFF3 Overlap
# ============

immutable OverlapIterator
    reader::Reader
    interval::Interval
end

function Base.eltype(::Type{OverlapIterator})
    return Record
end

function Base.iteratorsize(::Type{OverlapIterator})
    return Base.SizeUnknown()
end

function GenomicFeatures.eachoverlap(reader::Reader, interval::Interval)
    if isnull(reader.index)
        throw(ArgumentError("index is null"))
    end
    return OverlapIterator(reader, interval)
end

type OverlapIteratorState
    chunks::Vector{GenomicFeatures.Indexes.Chunk}
    chunkid::Int
    record::Record
    done::Bool
end

function Base.start(iter::OverlapIterator)
    @assert !isnull(iter.reader.index)
    return OverlapIteratorState(
        GenomicFeatures.Indexes.overlapchunks(get(iter.reader.index), iter.interval),
        0,
        Record(),
        false)
end

function Base.done(iter::OverlapIterator, state)
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

function Base.next(iter::OverlapIterator, state)
    return copy(state.record), state
end

function icmp(record, interval)
    c = cmp(seqid(record), interval.seqname)
    if c < 0 || (c == 0 && seqend(record) < interval.first)
        return -1
    elseif c > 0 || (c == 0 && seqstart(record) > interval.last)
        return +1
    else
        return 0
    end
end
