# Interval
# ========
#
# Base interval types and utilities.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

# Note, just to be clear: this shadows IntervalTrees.Interval
"""
    struct Interval{T} <: IntervalTrees.AbstractInterval{Int64}

The first three fields (`seqname`, `first`, and `last`) are mandatory arguments when constructing the [`Interval`](@ref Interval) object.

# Fields
- `seqname::String`: the sequence name associated with the interval.
- `first::Int64`: the leftmost position.
- `last::Int64`: the rightmost position.
- `strand::Strand`: the [`strand`](@ref Strand).
- `metadata::T`
"""
struct Interval{T} <: IntervalTrees.AbstractInterval{Int64}
    seqname::String
    first::Int64
    last::Int64
    strand::Strand
    metadata::T
end

function Interval(seqname::AbstractString, first::Integer, last::Integer, strand::Union{Strand,Char}=STRAND_BOTH, metadata::T=nothing) where T
    return Interval{T}(seqname, first, last, strand, metadata)
end

function Interval(seqname::AbstractString, range::UnitRange{R}, strand::Union{Strand,Char}=STRAND_BOTH, metadata::T=nothing) where {T,R<:Integer}
    return Interval{T}(seqname, first(range), last(range), strand, metadata)
end


"""
    Interval{T}(data)

The returned data is converted to Interval{T} if there is an implemented [`Base.convert`](https://docs.julialang.org/en/v1/base/base/#Base.convert) function for the type of data.
This method provides a useful hook for converting custom types to Interval{T}.
"""
function Interval{T}(data) :: Interval{T} where T
    return data #Note: the returned data is converted to Interval{T}.
end

function BioGenerics.seqname(i::Interval)
    return i.seqname
end

function BioGenerics.metadata(i::Interval)
    return i.metadata
end

function strand(i::Interval)
    return i.strand
end

"""
    leftposition(i::Interval)

Return the leftmost position of `i`.
"""
function BioGenerics.leftposition(i::Interval)
    return i.first
end

"""
    rightposition(i::Interval)

Return the rightmost position of `i`.
"""
function BioGenerics.rightposition(i::Interval)
    return i.last
end

IntervalTrees.first(i::Interval) = leftposition(i)
IntervalTrees.last(i::Interval) = rightposition(i)

function Base.isless(a::Interval, b::Interval, seqname_isless::Function=isless)
    if seqname(a) != seqname(b)
        return seqname_isless(seqname(a), seqname(b))::Bool
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

Intervals are considered well ordered if seqname(a) <= seqname(b) and leftposition(a) <= leftposition(b).
"""
function isordered(a::Interval, b::Interval, seqname_isless::Function=isless)
    if seqname(a) != seqname(b)
        return seqname_isless(seqname(a), seqname(b))::Bool
    end

    if leftposition(a) != leftposition(b)
        return leftposition(a) < leftposition(b)
    end

    return true
end

"""
Return true if interval `a` entirely precedes `b`.
"""
function precedes(a::Interval, b::Interval, seqname_isless::Function=isless)
    return (rightposition(a) < leftposition(b) && seqname(a) == seqname(b)) || seqname_isless(seqname(a), seqname(b))::Bool
end

function Base.:(==)(a::Interval, b::Interval)
    return seqname(a)       == seqname(b) &&
           leftposition(a)  == leftposition(b) &&
           rightposition(a) == rightposition(b) &&
           strand(a)        == strand(b) &&
           metadata(a)      == metadata(b)
end

"Return true if interval `a` overlaps interval `b`, with no consideration to strand"
function BioGenerics.isoverlapping(a::Interval, b::Interval)
    return leftposition(a) <= rightposition(b) &&
           leftposition(b) <= rightposition(a) &&
           seqname(a)      == seqname(b)
end

function Base.show(io::IO, i::Interval)
    if get(io, :compact, false)
        print(io, seqname(i), ":", leftposition(i), "-", rightposition(i), "  ", strand(i), "  ", metadata(i) === nothing ? "nothing" : metadata(i))
    else
        println(io, summary(i), ':')
        println(io, "  sequence name: ", seqname(i))
        println(io, "  leftmost position: ", leftposition(i))
        println(io, "  rightmost position: ", rightposition(i))
        println(io, "  strand: ", strand(i))
          print(io, "  metadata: ", metadata(i) === nothing ? "nothing" : metadata(i))
    end
end

function metadatatype(::Type{T}) where T
    return _metadatatype(eltype(T))
end

function metadatatype(x::Any)
    return metadatatype(typeof(x))
end

function _metadatatype(::Type{Interval{T}}) where T
    return T
end

function _metadatatype(::Type{T}) where T
    return T
end
