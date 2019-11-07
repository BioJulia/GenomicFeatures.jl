# GenomicInterval
# ========
#
# Base interval types and utilities.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

abstract type AbstractGenomicInterval{T} <: IntervalTrees.AbstractInterval{Int64} end

"A genomic interval specifies interval with some associated metadata."
struct GenomicInterval{T} <: AbstractGenomicInterval{T}
    seqname::String
    first::Int64
    last::Int64
    strand::Strand
    metadata::T
end

function GenomicInterval(seqname::AbstractString, first::Integer, last::Integer, strand::Union{Strand,Char}=STRAND_BOTH, metadata=nothing)
    return GenomicInterval{typeof(metadata)}(seqname, first, last, strand, metadata)
end

function GenomicInterval(seqname::AbstractString, range::UnitRange{T}, strand::Union{Strand,Char}=STRAND_BOTH, metadata=nothing) where T<:Integer
    return GenomicInterval{typeof(metadata)}(seqname, first(range), last(range), strand, metadata)
end

function BioGenerics.seqname(i::AbstractGenomicInterval)
    return i.seqname
end

function BioGenerics.metadata(i::AbstractGenomicInterval)
    return i.metadata
end

function strand(i::GenomicInterval)
    return i.strand
end

"""
    leftposition(i::AbstractGenomicInterval)

Return the leftmost position of `i`.
"""
function BioGenerics.leftposition(i::AbstractGenomicInterval)
    return i.first
end

"""
    rightposition(i::AbstractGenomicInterval)

Return the rightmost position of `i`.
"""
function BioGenerics.rightposition(i::AbstractGenomicInterval)
    return i.last
end

IntervalTrees.first(i::AbstractGenomicInterval) = leftposition(i)
IntervalTrees.last(i::AbstractGenomicInterval) = rightposition(i)

function Base.isless(a::AbstractGenomicInterval{T}, b::AbstractGenomicInterval{T}, seqname_isless::Function=isless) where T

    a_seqname = seqname(a)
    b_seqname = seqname(b)

    if a_seqname != b_seqname
        return seqname_isless(a_seqname, b_seqname)::Bool
    end

    a_leftposition = leftposition(a)
    b_leftposition = leftposition(b)

    if a_leftposition != b_leftposition
        return a_leftposition < b_leftposition
    end

    a_rightposition = rightposition(a)
    b_rightposition = rightposition(b)

    if a_rightposition != b_rightposition
        return a_rightposition < b_rightposition
    end

    return false
end

function Base.isless(a::GenomicInterval{T}, b::GenomicInterval{T}, seqname_isless::Function=isless) where T
    if a.seqname != b.seqname
        return seqname_isless(a.seqname, b.seqname)::Bool
    elseif a.first != b.first
        return a.first < b.first
    elseif a.last != b.last
        return a.last < b.last
    elseif a.strand != b.strand
        return a.strand < b.strand
    else
        return false
    end
end

"""
Check if two intervals are well ordered.

AbstractGenomicIntervals are considered well ordered if seqname(a) <= seqname(b)nd and leftposition(a) <= leftposition(b).
"""
function isordered(a::AbstractGenomicInterval{T}, b::AbstractGenomicInterval{T}, seqname_isless::Function=isless) where T

    a_seqname = seqname(a)
    b_seqname = seqname(b)

    if a_seqname != b_seqname
        return seqname_isless(a_seqname, b_seqname)::Bool
    end

    a_leftposition = leftposition(a)
    b_leftposition = leftposition(b)

    if a_leftposition != b_leftposition
        return a_leftposition < b_leftposition
    end

    return true
end

"""
Return true if interval `a` entirely precedes `b`.
"""
function precedes(a::AbstractGenomicInterval{T}, b::AbstractGenomicInterval{T}, seqname_isless::Function=isless) where T
    return (rightposition(a) < leftposition(b) && seqname(a) == seqname(b)) || seqname_isless(seqname(a), seqname(b))::Bool
end

function Base.:(==)(a::AbstractGenomicInterval{T}, b::AbstractGenomicInterval{T}) where T
    return seqname(a) == seqname(b) &&
           leftposition(a) == leftposition(b) &&
           rightposition(a) == rightposition(b)
end

function Base.:(==)(a::GenomicInterval{T}, b::GenomicInterval{T}) where T
    return a.seqname  == b.seqname &&
           a.first    == b.first &&
           a.last     == b.last &&
           a.strand   == b.strand &&
           a.metadata == b.metadata
end

"Return true if interval `a` overlaps interval `b`, with no consideration to strand"
function BioGenerics.isoverlapping(a::AbstractGenomicInterval{S}, b::AbstractGenomicInterval{T}) where {S, T}
    return leftposition(a) <= rightposition(b) && leftposition(b) <= rightposition(a) && seqname(a) == seqname(b)
end

function Base.show(io::IO, i::AbstractGenomicInterval)
    if get(io, :compact, false)
        print(io, seqname(i), ":", leftposition(i), "-", rightposition(i))
    else
        println(io, summary(i), ':')
        println(io, "  sequence name: ", seqname(i))
        println(io, "  leftmost position: ", leftposition(i))
          print(io, "  rightmost position: ", rightposition(i))
    end
end

function Base.show(io::IO, i::GenomicInterval)
    if get(io, :compact, false)
        print(io, i.seqname, ":", i.first, "-", i.last, "  ", i.strand, "  ", i.metadata === nothing ? "nothing" : i.metadata)
    else
        println(io, summary(i), ':')
        println(io, "  sequence name: ", i.seqname)
        println(io, "  leftmost position: ", i.first)
        println(io, "  rightmost position: ", i.last)
        println(io, "  strand: ", i.strand)
          print(io, "  metadata: ", i.metadata === nothing ? "nothing" : i.metadata)
    end
end

function metadatatype(::Type{T}) where T
    return _metadatatype(eltype(T))
end

function metadatatype(x::Any)
    return metadatatype(typeof(x))
end

function _metadatatype(::Type{A}) where {T, A <: AbstractGenomicInterval{T}}
    return T
end

function _metadatatype(::Type{T}) where T
    return T
end
