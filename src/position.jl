abstract type AbstractGenomicPosition{T} <: GenomicFeatures.AbstractGenomicInterval{T} end

"A genomic position is distinguished from an interval and only specifies a position with some associated metadata."
mutable struct GenomicPosition{T} <: AbstractGenomicPosition{T}
    groupname::String
    pos::Int64
    metadata::T
end

function GenomicPosition(groupname::String, pos::Integer, metadata::T = nothing) where {T}
    return GenomicPosition{T}(groupname, pos, metadata)
end

function GenomicPosition(groupname::String, range::UnitRange{<:Integer}, metadata::T = nothing) where {T}
    length(range) == 1 || error("Position must have a length of 1.")
    return GenomicPosition(groupname, first(range), metadata)
end

function left_position(site::AbstractGenomicPosition)
    return site.pos
end

function right_position(site::AbstractGenomicPosition)
    return site.pos
end

"""
    position(p::AbstractGenomicPosition)

Return the position of `p`.
"""
function Base.position(site::AbstractGenomicPosition)
    return left_position(site)
end

function Base.show(io::IO, p::AbstractGenomicPosition)
    return if get(io, :compact, false)
        print(io, group_name(p), ":", position(p))
    else
        println(io, summary(p), ':')
        println(io, "  group name: ", group_name(p))
        print(io, "  position: ", position(p))
    end
end

function Base.show(io::IO, p::GenomicPosition)
    return if get(io, :compact, false)
        print(io, group_name(p), ":", position(p), "  ", metadata(p) === nothing ? "nothing" : metadata(p))
    else
        println(io, summary(p), ':')
        println(io, "  group name: ", group_name(p))
        println(io, "  position: ", position(p))
        print(io, "  metadata: ", metadata(p) === nothing ? "nothing" : metadata(p))
    end
end

function Base.:(==)(a::GenomicPosition, b::GenomicPosition)
    return group_name(a) == group_name(b) &&
        position(a) == position(b) &&
        metadata(a) == metadata(b)
end
