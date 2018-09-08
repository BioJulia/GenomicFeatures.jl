# Queue
# =====

mutable struct Queue{T}
    data::Vector{T}
    offset::Int
    first::Int
    last::Int
end

function Queue{T}(bufsize::Integer=2^4) where {T}
    if bufsize ≤ 0
        throw(ArgumentError("buffer size must be positive"))
    elseif !ispow2(bufsize)
        throw(ArgumentError("buffer size must be a power of two"))
    end
    return Queue(Vector{T}(undef, bufsize), 0, 1, 0)
end

function Base.eltype(::Type{Queue{T}}) where T
    return T
end

function Base.length(queue::Queue)
    return queue.last - queue.first + 1
end

function Base.isempty(queue::Queue)
    return queue.first > queue.last
end

function Base.push!(queue::Queue{T}, elm::T) where T
    if length(queue.data) < length(queue) + 1
        index_first = dataindex(queue, queue.first)
        index_last = dataindex(queue, queue.last)
        index_end = lastindex(queue.data)
        # NOTE: resize factor must be a power of two
        resize!(queue.data, 2 * length(queue.data))
        @assert ispow2(length(queue.data))
        if !isempty(queue) && index_last < index_first
            # make the circular data linear
            copyto!(queue.data, index_end + 1, queue.data, 1, index_last)
        end
        copyto!(queue.data, 1, queue.data, index_first, length(queue))
        queue.offset = 0
    end
    queue.data[dataindex(queue, queue.last + 1)] = elm
    queue.last += 1
    return queue
end

function Base.popfirst!(queue::Queue)
    if isempty(queue)
        throw(ArgumentError("empty"))
    end
    elm = queue.data[dataindex(queue, queue.first)]
    queue.first += 1
    queue.offset += 1
    return elm
end

function Base.iterate(queue::Queue)
    return iterate(queue, queue.first)
end

function Base.iterate(queue::Queue, i)
    if !(queue.first <= i <= queue.last)
        return nothing
    else
        queue.data[dataindex(queue, i)], i + 1
    end
end

function dataindex(queue::Queue, i::Integer)
    return ((i - queue.first + queue.offset) & (length(queue.data) - 1)) + 1
end

function Base.firstindex(queue::Queue)
    return queue.first
end

function Base.lastindex(queue::Queue)
    return queue.last
end

function Base.getindex(queue::Queue, i::Integer)
    if !(queue.first <= i <= queue.last)
        throw(BoundsError())
    end

    return queue.data[dataindex(queue, i)]
end


