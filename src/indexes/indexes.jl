# Index
# =====
#
# Index types for genomic intervals.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

module Indexes

import BGZFStreams: BGZFStream, VirtualOffset
import GenomicFeatures: Interval

include("chunk.jl")
include("bgzfindex.jl")
include("tabix.jl")

end
