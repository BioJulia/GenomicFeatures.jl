# Tabix Overlap Iterator
# ======================

struct TabixOverlapIterator{T}
    reader::T
    interval::Interval
end

function Base.eltype(::Type{TabixOverlapIterator{T}}) where T
    return eltype(T)
end

function Base.IteratorSize(::Type{TabixOverlapIterator{T}}) where T
    return Base.SizeUnknown()
end

mutable struct TabixOverlapIteratorState{T}
    chunks::Vector{Indexes.Chunk}
    chunkid::Int
    done::Bool
    record::T
end

function Base.iterate(iter::TabixOverlapIterator)
    @assert iter.reader.index !== nothing
    # TODO: Use a method that resets the reading position.
    buffer = BioCore.IO.stream(iter.reader)
    iter.reader.state = BioCore.Ragel.State(1, BufferedStreams.BufferedInputStream(buffer.source))
    state = TabixOverlapIteratorState(
        Indexes.overlapchunks(iter.reader.index, iter.interval),
        0,
        false,
        eltype(iter)())

    return iterate(iter, state)
end

function done(iter::TabixOverlapIterator, state)
    buffer = BioCore.IO.stream(iter.reader)
    source = buffer.source
    if state.chunkid == 0
        if isempty(state.chunks)
            return true
        end
        state.chunkid += 1
        seek(source, state.chunks[state.chunkid].start)
    end
    while state.chunkid ≤ lastindex(state.chunks)
        chunk = state.chunks[state.chunkid]
        # The `virtualoffset(source)` is not synchronized with the current
        # reading position because data are buffered in `buffer` for parsing
        # text. So we need to check not only `virtualoffset` but also
        # `nb_available`, which returns the current buffered data size.
        while bytesavailable(buffer) > 0 || BGZFStreams.virtualoffset(source) < chunk.stop
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
        if state.chunkid ≤ lastindex(state.chunks)
            seek(source, state.chunks[state.chunkid].start)
        end
    end
    # no more overlapping records
    return true
end

function Base.iterate(iter::TabixOverlapIterator, state)
    if done(iter, state)
        return nothing
    else
        return copy(state.record), state
    end
end

function icmp(record, interval)
    c = cmp(BioCore.seqname(record), interval.seqname)
    if c < 0 || (c == 0 && BioCore.rightposition(record) < interval.first)
        return -1
    elseif c > 0 || (c == 0 && BioCore.leftposition(record) > interval.last)
        return +1
    else
        return 0
    end
end
