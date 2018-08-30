__precompile__()

module GenomicFeatures

export
    Strand,
    STRAND_NA,
    STRAND_POS,
    STRAND_NEG,
    STRAND_BOTH,

    Interval,
    seqname,
    leftposition,
    rightposition,
    strand,
    metadata,
    isoverlapping,

    IntervalCollection,
    eachoverlap,
    coverage,

    BED,
    GFF3,
    BigWig,
    BigBed,
    isfilled,
    hasseqname,
    hasleftposition,
    hasrightposition

import BGZFStreams
import BioCore: BioCore, seqname, leftposition, rightposition, isoverlapping, isfilled, hasseqname, hasleftposition, hasrightposition, metadata
import DataStructures
import IntervalTrees
import Base.@propagate_inbounds

include("strand.jl")
include("interval.jl")
include("intervalcollection.jl")
include("queue.jl")
include("overlap.jl")
include("coverage.jl")
include("indexes/indexes.jl")
include("bed/bed.jl")
include("gff3/gff3.jl")
include("bbi/bbi.jl")
include("bigwig/bigwig.jl")
include("bigbed/bigbed.jl")

end # module
