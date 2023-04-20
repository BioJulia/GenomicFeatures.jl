# GenomicInterval
# ========
#
# Base interval types and utilities.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

"""
    GenomicInterval{T} <: AbstractGenomicInterval{T}

A genomic interval specifies interval with some associated metadata.
The first three fields (`groupname`, `first`, and `last`) are mandatory arguments when constructing the [`Interval`](@ref Interval) object.

# Fields
- `groupname::String`: the group name associated with the interval.
- `first::Int64`: the leftmost position.
- `last::Int64`: the rightmost position.
- `strand::Strand`: the [`strand`](@ref Strand).
- `metadata::T`
"""
struct GenomicInterval{T} <: AbstractGenomicInterval{T}
    groupname::String
    first::Int64
    last::Int64
    strand::Strand
    metadata::T
end

function GenomicInterval(groupname::AbstractString, first::Integer, last::Integer, strand::Union{Strand,Char}=STRAND_BOTH, metadata::T=nothing) where T
    return GenomicInterval{T}(groupname, first, last, strand, metadata)
end

function GenomicInterval(groupname::AbstractString, range::UnitRange{R}, strand::Union{Strand,Char}=STRAND_BOTH, metadata::T=nothing) where {T,R<:Integer}
    return GenomicInterval{T}(groupname, first(range), last(range), strand, metadata)
end


"""
    GenomicInterval{T}(data)

The returned data is converted to GenomicInterval{T} if there is an implemented [`Base.convert`](https://docs.julialang.org/en/v1/base/base/#Base.convert) function for the type of data.
This method provides a useful hook for converting custom types to GenomicInterval{T}.
"""
function GenomicInterval{T}(data) :: GenomicInterval{T} where T
    return data #Note: the returned data is converted to GenomicInterval{T}.
end

function BioGenerics.groupname(i::AbstractGenomicInterval)
    return i.groupname
end

function BioGenerics.seqname(i::AbstractGenomicInterval)
    return groupname(i)
end

function BioGenerics.metadata(i::AbstractGenomicInterval)
    return i.metadata
end

function BioGenerics.metadata(i::AbstractGenomicInterval{Nothing})
    return nothing
end

function strand(i::AbstractGenomicInterval)
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

function Base.isless(a::AbstractGenomicInterval, b::AbstractGenomicInterval, groupname_isless::Function=isless)
    if groupname(a) != groupname(b)
        return groupname_isless(groupname(a), groupname(b))::Bool
    end

    if leftposition(a) != leftposition(b)
        return leftposition(a) < leftposition(b)
    end

    if rightposition(a) != rightposition(b)
        return rightposition(a) < rightposition(b)
    end

    return false
end

function Base.isless(a::GenomicInterval, b::GenomicInterval, groupname_isless::Function=isless)
    if groupname(a) != groupname(b)
        return groupname_isless(groupname(a), groupname(b))::Bool
    end

    if leftposition(a) != leftposition(b)
        return leftposition(a) < leftposition(b)
    end

    if rightposition(a) != rightposition(b)
        return rightposition(a) < rightposition(b)
    end

    if strand(a) != strand(b)
        return strand(a) < strand(b)
    end

    return false
end

"""
Check if two intervals are well ordered.

`AbstractGenomicInterval` are considered well ordered if groupname(a) <= groupname(b) and leftposition(a) <= leftposition(b).
"""
function isordered(a::AbstractGenomicInterval, b::AbstractGenomicInterval, groupname_isless::Function=isless)
    if groupname(a) != groupname(b)
        return groupname_isless(groupname(a), groupname(b))::Bool
    end

    if leftposition(a) != leftposition(b)
        return leftposition(a) < leftposition(b)
    end

    return true
end

"""
Return true if interval `a` entirely precedes `b`.
"""
function precedes(a::AbstractGenomicInterval, b::AbstractGenomicInterval, groupname_isless::Function=isless)

    if groupname(a) != groupname(b)
        return groupname_isless(groupname(a), groupname(b))::Bool
    end

    return rightposition(a) < leftposition(b)
end

function Base.:(==)(a::AbstractGenomicInterval, b::AbstractGenomicInterval)
    return groupname(a)     == groupname(b) &&
           leftposition(a)  == leftposition(b) &&
           rightposition(a) == rightposition(b) &&
           metadata(a)      == metadata(b)
end

function Base.:(==)(a::GenomicInterval, b::GenomicInterval)
    return groupname(a)     == groupname(b) &&
           leftposition(a)  == leftposition(b) &&
           rightposition(a) == rightposition(b) &&
           strand(a)        == strand(b) &&
           metadata(a)      == metadata(b)
end

"Return true if interval `a` overlaps interval `b`, with no consideration to strand"
function BioGenerics.isoverlapping(a::AbstractGenomicInterval, b::AbstractGenomicInterval)
    return leftposition(a) <= rightposition(b) &&
           leftposition(b) <= rightposition(a) &&
           groupname(a)    == groupname(b)
end

function Base.show(io::IO, i::AbstractGenomicInterval)
    if get(io, :compact, false)
        print(io, groupname(i), ":", leftposition(i), "-", rightposition(i), "  ", metadata(i) === nothing ? "nothing" : metadata(i))
    else
        println(io, summary(i), ':')
        println(io, "  group name: ", groupname(i))
        println(io, "  leftmost position: ", leftposition(i))
        println(io, "  rightmost position: ", rightposition(i))
          print(io, "  metadata: ", metadata(i) === nothing ? "nothing" : metadata(i))
    end
end

function Base.show(io::IO, i::GenomicInterval)
    if get(io, :compact, false)
        print(io, groupname(i), ":", leftposition(i), "-", rightposition(i), "  ", strand(i), "  ", metadata(i) === nothing ? "nothing" : metadata(i))
    else
        println(io, summary(i), ':')
        println(io, "  group name: ", groupname(i))
        println(io, "  leftmost position: ", leftposition(i))
        println(io, "  rightmost position: ", rightposition(i))
        println(io, "  strand: ", strand(i))
          print(io, "  metadata: ", metadata(i) === nothing ? "nothing" : metadata(i))
    end
end
