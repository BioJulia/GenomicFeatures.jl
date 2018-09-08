# GFF3 File Format
# ================

module GFF3

import Automa
import Automa.RegExp: @re_str
import BGZFStreams
import BioCore.Exceptions: missingerror
import BioSequences
import BufferedStreams
import GenomicFeatures: GenomicFeatures, Interval, IntervalCollection
import URIParser
using BioCore

include("record.jl")
include("reader.jl")
include("writer.jl")

end
