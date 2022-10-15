__precompile__()

module GenomicFeatures

export
    Strand,
    STRAND_NA,
    STRAND_POS,
    STRAND_NEG,
    STRAND_BOTH,

    GenomicInterval,
	groupname,
    seqname,
    leftposition,
    rightposition,
    strand,
    metadata,
    isoverlapping,

    GenomicPosition,

    GenomicIntervalCollection,
    eachoverlap,
    coverage,
	hasintersection,

    isfilled,
    hasseqname,
    hasleftposition,
    hasrightposition,

	span

import BioGenerics: BioGenerics, groupname, seqname, leftposition, rightposition, isoverlapping, isfilled, hasgroupname, hasseqname, hasleftposition, hasrightposition, metadata
import DataStructures
import IntervalTrees
import Base.@propagate_inbounds

abstract type AbstractGenomicInterval{T} <: IntervalTrees.AbstractInterval{Int64} end

abstract type AbstractGenomicCollection{I} end

const _IterableGenomicCollection{I} = Union{<:AbstractVector{I}, <:AbstractGenomicCollection{I}} where {I}
const IterableGenomicCollection{I} = Union{<:_IterableGenomicCollection{I}, <:Base.Generator{<:_IterableGenomicCollection{I}, <:Any}} where {I}

function Base.eltype(::Type{C}) where {I, C<:Base.Generator{<:_IterableGenomicCollection{I}, <:Any}}
    return I
end

include("strand.jl")
include("interval.jl")
include("position.jl")
include("intervalcollection.jl")
include("queue.jl")
include("overlap.jl")
include("coverage.jl")

"""
    span(interval::AbstractGenomicInterval)::Int

Get the span of `interval`.
"""
function span(interval::AbstractGenomicInterval)
	return length(leftposition(interval):rightposition(interval))
end

"""
    volume(interval::AbstractGenomicInterval)

Get the product of the `interval`'s span and metadata.
"""
function volume(interval::AbstractGenomicInterval)
	return span(interval) * GenomicFeatures.metadata(interval)
end

function intervaltype(::Type{I}) where {I<:AbstractGenomicInterval}
    return I
end

"""
Get the interval type.
Overwrite to suggest stream element conversions during iteration.
"""
function intervaltype(::Type{T}) where T
    return GenomicInterval{T}
end

function intervaltype(interval::T) where T
    return intervaltype(I)
end

"""
Get the interval type without metadata.
This query is useful when converting an intervals metadata.
"""
function baseintervaltype(::Type{I}) where {I<:GenomicInterval}
	return GenomicInterval
end

function baseintervaltype(::Type{I}) where {I<:GenomicPosition}
    return GenomicPosition
end

function baseintervaltype(::I) where {I<:GenomicFeatures.AbstractGenomicInterval} #TODO: consider loosening captured types.
	return baseintervaltype(I)
end

"""
Get the interval's metadata type.
"""
function intervalmetadatatype(interval::Type{I}) where {T, I<:GenomicFeatures.AbstractGenomicInterval{T}}
    return T
end

function intervalmetadatatype(::I) where {T, I<:GenomicFeatures.AbstractGenomicInterval{T}} #TODO: consider loosening captured types.
    return intervalmetadatatype(I)
end

end # module
