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
The first three fields (`groupname`, `first`, and `last`) are mandatory arguments when constructing the [`GenomicInterval`](@ref) object.

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

function GenomicInterval(groupname::AbstractString, first::Integer, last::Integer, strand::Union{Strand, Char} = STRAND_BOTH, metadata::T = nothing) where {T}
    return GenomicInterval{T}(groupname, first, last, strand, metadata)
end

function GenomicInterval(groupname::AbstractString, range::UnitRange{R}, strand::Union{Strand, Char} = STRAND_BOTH, metadata::T = nothing) where {T, R <: Integer}
    return GenomicInterval{T}(groupname, first(range), last(range), strand, metadata)
end


"""
    GenomicInterval{T}(data)

The returned data is converted to GenomicInterval{T} if there is an implemented [`Base.convert`](https://docs.julialang.org/en/v1/base/base/#Base.convert) function for the type of data.
This method provides a useful hook for converting custom types to GenomicInterval{T}.
"""
function GenomicInterval{T}(data)::GenomicInterval{T} where {T}
    return data #Note: the returned data is converted to GenomicInterval{T}.
end

function group_name(i::AbstractGenomicInterval)
    return i.groupname
end

function BioGenerics.seqname(i::AbstractGenomicInterval)
    return group_name(i)
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
    left_position(i::AbstractGenomicInterval)

Return the leftmost position of `i`.
"""
function left_position(i::AbstractGenomicInterval)
    return i.first
end

"""
    right_position(i::AbstractGenomicInterval)

Return the rightmost position of `i`.
"""
function right_position(i::AbstractGenomicInterval)
    return i.last
end

IntervalTrees.first(i::AbstractGenomicInterval) = left_position(i)
IntervalTrees.last(i::AbstractGenomicInterval) = right_position(i)

function Base.isless(a::AbstractGenomicInterval, b::AbstractGenomicInterval, groupname_isless::Function = isless)
    if group_name(a) != group_name(b)
        return groupname_isless(group_name(a), group_name(b))::Bool
    end

    if left_position(a) != left_position(b)
        return left_position(a) < left_position(b)
    end

    if right_position(a) != right_position(b)
        return right_position(a) < right_position(b)
    end

    return false
end

function Base.isless(a::GenomicInterval, b::GenomicInterval, groupname_isless::Function = isless)
    if group_name(a) != group_name(b)
        return groupname_isless(group_name(a), group_name(b))::Bool
    end

    if left_position(a) != left_position(b)
        return left_position(a) < left_position(b)
    end

    if right_position(a) != right_position(b)
        return right_position(a) < right_position(b)
    end

    if strand(a) != strand(b)
        return strand(a) < strand(b)
    end

    return false
end

"""
Check if two intervals are well ordered.

`AbstractGenomicInterval` are considered well ordered if group_name(a) <= group_name(b) and left_position(a) <= left_position(b).
"""
function is_ordered(a::AbstractGenomicInterval, b::AbstractGenomicInterval, groupname_isless::Function = isless)
    if group_name(a) != group_name(b)
        return groupname_isless(group_name(a), group_name(b))::Bool
    end

    if left_position(a) != left_position(b)
        return left_position(a) < left_position(b)
    end

    return true
end

"""
Return true if interval `a` entirely precedes `b`.
"""
function precedes(a::AbstractGenomicInterval, b::AbstractGenomicInterval, groupname_isless::Function = isless)

    if group_name(a) != group_name(b)
        return groupname_isless(group_name(a), group_name(b))::Bool
    end

    return right_position(a) < left_position(b)
end

function Base.:(==)(a::AbstractGenomicInterval, b::AbstractGenomicInterval)
    return group_name(a) == group_name(b) &&
        left_position(a) == left_position(b) &&
        right_position(a) == right_position(b) &&
        metadata(a) == metadata(b)
end

function Base.:(==)(a::GenomicInterval, b::GenomicInterval)
    return group_name(a) == group_name(b) &&
        left_position(a) == left_position(b) &&
        right_position(a) == right_position(b) &&
        strand(a) == strand(b) &&
        metadata(a) == metadata(b)
end

"Return true if interval `a` overlaps interval `b`, with no consideration to strand"
function is_overlapping(a::AbstractGenomicInterval, b::AbstractGenomicInterval)
    return left_position(a) <= right_position(b) &&
        left_position(b) <= right_position(a) &&
        group_name(a) == group_name(b)
end

function Base.show(io::IO, i::AbstractGenomicInterval)
    return if get(io, :compact, false)
        print(io, group_name(i), ":", left_position(i), "-", right_position(i), "  ", metadata(i) === nothing ? "nothing" : metadata(i))
    else
        println(io, summary(i), ':')
        println(io, "  group name: ", group_name(i))
        println(io, "  leftmost position: ", left_position(i))
        println(io, "  rightmost position: ", right_position(i))
        print(io, "  metadata: ", metadata(i) === nothing ? "nothing" : metadata(i))
    end
end

function Base.show(io::IO, i::GenomicInterval)
    return if get(io, :compact, false)
        print(io, group_name(i), ":", left_position(i), "-", right_position(i), "  ", strand(i), "  ", metadata(i) === nothing ? "nothing" : metadata(i))
    else
        println(io, summary(i), ':')
        println(io, "  group name: ", group_name(i))
        println(io, "  leftmost position: ", left_position(i))
        println(io, "  rightmost position: ", right_position(i))
        println(io, "  strand: ", strand(i))
        print(io, "  metadata: ", metadata(i) === nothing ? "nothing" : metadata(i))
    end
end
