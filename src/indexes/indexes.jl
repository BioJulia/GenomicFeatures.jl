# Index
# =====
#
# Index types for genomic intervals.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

module Indexes

import BGZFStreams
import BioCore
import BufferedStreams
import GenomicFeatures: Interval

include("chunk.jl")
include("bgzfindex.jl")
include("tabix.jl")
include("overlap.jl")

end
