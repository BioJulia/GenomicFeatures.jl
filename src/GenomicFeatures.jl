__precompile__()

module GenomicFeatures

    export
        Strand,
        STRAND_NA,
        STRAND_POS,
        STRAND_NEG,
        STRAND_BOTH,

        GenomicInterval,
        group_name,
        left_position,
        right_position,
        strand,
        metadata,
        is_overlapping,

        GenomicPosition,

        GenomicIntervalCollection,
        each_overlap,
        coverage,
        has_intersection,

        is_filled,
        has_group_name,
        has_left_position,
        has_right_position,

        span

    # The squashed names of GenomicFeatures 3.0 are deprecated bindings, which `deprecated.jl` defines and exports.
    import BioGenerics
    import BioGenerics: metadata
    import DataStructures
    import IntervalTrees
    import Base.@propagate_inbounds

    # The snake case names of the functions that BioGenerics owns.
    include("names.jl")

    abstract type AbstractGenomicInterval{T} <: IntervalTrees.AbstractInterval{Int64} end

    abstract type AbstractGenomicCollection{I} end

    const _IterableGenomicCollection{I} = Union{<:AbstractVector{I}, <:AbstractGenomicCollection{I}} where {I}
    const IterableGenomicCollection{I} = Union{<:_IterableGenomicCollection{I}, <:Base.Generator{<:_IterableGenomicCollection{I}, <:Any}} where {I}

    function Base.eltype(::Type{C}) where {I, C <: Base.Generator{<:_IterableGenomicCollection{I}, <:Any}}
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
        return length(left_position(interval):right_position(interval))
    end

    """
        volume(interval::AbstractGenomicInterval)

    Get the product of the `interval`'s span and metadata.
    """
    function volume(interval::AbstractGenomicInterval)
        return span(interval) * GenomicFeatures.metadata(interval)
    end

    """
    Get the interval type.
    Overwrite to suggest stream element conversions during iteration.
    """
    function interval_type(::Type{I}) where {I <: AbstractGenomicInterval}
        return I
    end

    function interval_type(interval::T) where {T}
        return interval_type(T)
    end

    """
    Get the interval type without metadata.
    This query is useful when converting an interval's metadata or when implementing promotion rules.
    """
    function base_interval_type(::Type{I}) where {I <: GenomicFeatures.AbstractGenomicInterval} #TODO: consider loosening captured types to IntervalTrees.AbstractInterval.
        tn = Base.typename(I)
        return getproperty(tn.module, tn.name)
    end

    function base_interval_type(::I) where {I <: GenomicFeatures.AbstractGenomicInterval} #TODO: consider loosening captured types to IntervalTrees.AbstractInterval.
        return base_interval_type(I)
    end

    """
    Get the interval's metadata type.
    """
    function metadata_type(interval::Type{I}) where {T, I <: GenomicFeatures.AbstractGenomicInterval{T}}
        return T
    end

    function metadata_type(::I) where {T, I <: GenomicFeatures.AbstractGenomicInterval{T}} #TODO: consider loosening captured types to IntervalTrees.AbstractInterval.
        return metadata_type(I)
    end

    # The deprecated squashed names of GenomicFeatures 3.0.
    include("deprecated.jl")

end # module
