import Base: @deprecate, depwarn

export IntervalCollection

# Interval
@deprecate Interval(seqname::AbstractString, first::Integer, last::Integer, strand::Union{Strand,Char}=STRAND_BOTH, metadata=nothing) GenomicInterval(seqname, first, last, strand, metadata)
@deprecate Interval(seqname::AbstractString, range::UnitRange{T}, strand::Union{Strand,Char}=STRAND_BOTH, metadata=nothing) where T<:Integer GenomicInterval(seqname, range, strand, metadata)

# IntervalCollection
struct IntervalCollection{T} # Work around for parametric types.
    function IntervalCollection{T}() where T
        depwarn(
            string(
                "The empty constructor IntervalCollection{$T}() is deprecated. ",
                "Defaulting to GenomicIntervalCollection{$T}()."
            ),
            :IntervalCollection #TODO: find a way to include parametric type.
        )

        return GenomicIntervalCollection{T}()
    end

    function IntervalCollection{T}(intervals::AbstractVector, sort::Bool=false) where T
        depwarn(
            string(
                "The constructor IntervalCollection{$T}(intervals, sort) is deprecated. ",
                "Defaulting to GenomicIntervalCollection{$T}(intervals, sort)."
            ),
            :IntervalCollection #TODO: find a way to include parametric type.
        )

        return GenomicIntervalCollection{T}(intervals, sort)
    end
end

@deprecate IntervalCollection GenomicIntervalCollection
# @deprecate IntervalCollection(intervals) GenomicIntervalCollection(intervals)
# @deprecate IntervalCollection{T} GenomicIntervalCollection{T} # Doesn't work with parametric types.
# @deprecate IntervalCollection{T}() where T GenomicIntervalCollection{T}() where T # Doesn't work with parametric types.
