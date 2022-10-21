__precompile__()

module GenomicFeatures

export
    Strand,
    STRAND_NA,
    STRAND_POS,
    STRAND_NEG,
    STRAND_BOTH,

    GenomicInterval,
    seqname,
    leftposition,
    rightposition,
    strand,
    metadata,
    isoverlapping,

    GenomicIntervalCollection,
    eachoverlap,
    coverage,
	hasintersection,

    isfilled,
    hasseqname,
    hasleftposition,
    hasrightposition,

	span

import BioGenerics: BioGenerics, seqname, leftposition, rightposition, isoverlapping, isfilled, hasseqname, hasleftposition, hasrightposition, metadata
import DataStructures
import IntervalTrees
import Base.@propagate_inbounds

include("strand.jl")
include("interval.jl")
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

end # module
