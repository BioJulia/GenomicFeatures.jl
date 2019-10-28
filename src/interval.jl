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

function BioGenerics.seqname(i::GenomicInterval)
    return i.seqname
end

function BioGenerics.metadata(i::GenomicInterval)
    return i.metadata
end

function strand(i::GenomicInterval)
    return i.strand
end

"""
    leftposition(i::GenomicInterval)

Return the leftmost position of `i`.
"""
function BioGenerics.leftposition(i::GenomicInterval)
    return i.first
end

"""
    rightposition(i::GenomicInterval)

Return the rightmost position of `i`.
"""
function BioGenerics.rightposition(i::GenomicInterval)
    return i.last
end

IntervalTrees.first(i::GenomicInterval) = i.first
IntervalTrees.last(i::GenomicInterval) = i.last

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

`GenomicIntervals` are considered well ordered if a.seqname <= b.seqnamend and a.first <= b.first.
"""
function isordered(a::GenomicInterval{T}, b::GenomicInterval{T}, seqname_isless::Function=isless) where T
    if a.seqname != b.seqname
        return seqname_isless(a.seqname, b.seqname)::Bool
    elseif a.first != b.first
        return a.first < b.first
    else
        return true
    end
end

"""
Return true if interval `a` entirely precedes `b`.
"""
function precedes(a::GenomicInterval{T}, b::GenomicInterval{T}, seqname_isless::Function=isless) where T
    return (a.last < b.first && a.seqname == b.seqname) || seqname_isless(a.seqname, b.seqname)::Bool
end

function Base.:(==)(a::GenomicInterval{T}, b::GenomicInterval{T}) where T
    return a.seqname  == b.seqname &&
           a.first    == b.first &&
           a.last     == b.last &&
           a.strand   == b.strand &&
           a.metadata == b.metadata
end

"Return true if interval `a` overlaps interval `b`, with no consideration to strand"
function BioGenerics.isoverlapping(a::GenomicInterval{S}, b::GenomicInterval{T}) where {S, T}
    return a.first <= b.last && b.first <= a.last && a.seqname == b.seqname
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

function _metadatatype(::Type{GenomicInterval{T}}) where T
    return T
end

function _metadatatype(::Type{T}) where T
    return T
end
