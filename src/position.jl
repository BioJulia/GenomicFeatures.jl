abstract type AbstractGenomicPosition{T} <: GenomicFeatures.AbstractGenomicInterval{T} end

"A genomic position is distinguished from an interval and only specifies a position with some associated metadata."
mutable struct GenomicPosition{T} <: AbstractGenomicPosition{T}
    seqname::String
    pos::Int64
    metadata::T
end

function GenomicPosition(seqname::String, pos::Integer, metadata::T = nothing) where T
    return GenomicPosition{T}(seqname, pos, metadata)
end

function GenomicPosition(seqname::String, range::UnitRange{<:Integer}, metadata::T = nothing) where T
    length(range) == 1 || error("Position must have a length of 1.")
    return GenomicPosition(seqname, first(range), metadata)
end

function BioGenerics.leftposition(site::AbstractGenomicPosition)
    return site.pos
end

function BioGenerics.rightposition(site::AbstractGenomicPosition)
    return site.pos
end

"""
    position(p::AbstractGenomicPosition)

Return the position of `p`.
"""
function Base.position(site::AbstractGenomicPosition)
    return leftposition(site)
end

function Base.show(io::IO, p::AbstractGenomicPosition)
    if get(io, :compact, false)
        print(io, seqname(p), ":", position(p))
    else
        println(io, summary(p), ':')
        println(io, "  sequence name: ", seqname(p))
          print(io, "  position: ", position(p))
    end
end

function Base.show(io::IO, p::GenomicPosition)
    if get(io, :compact, false)
        print(io, seqname(p), ":", position(p), "  ", metadata(p) === nothing ? "nothing" : metadata(p))
    else
        println(io, summary(p), ':')
        println(io, "  sequence name: ", seqname(p))
        println(io, "  position: ", position(p))
          print(io, "  metadata: ", metadata(p) === nothing ? "nothing" : metadata(p))
    end
end

function Base.:(==)(a::GenomicPosition, b::GenomicPosition) where T
    return seqname(a)  == seqname(b) &&
           position(a) == position(b) &&
           metadata(a) == metadata(b)
end
