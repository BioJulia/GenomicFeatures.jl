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

    GenomicPosition,

    GenomicIntervalCollection,
    eachoverlap,
    coverage,

    isfilled,
    hasseqname,
    hasleftposition,
    hasrightposition

import BioGenerics: BioGenerics, seqname, leftposition, rightposition, isoverlapping, isfilled, hasseqname, hasleftposition, hasrightposition, metadata
import DataStructures
import IntervalTrees
import Base.@propagate_inbounds

include("strand.jl")
include("interval.jl")
include("position.jl")
include("intervalcollection.jl")
include("queue.jl")
include("overlap.jl")
include("coverage.jl")

end # module
