# GFF3 File Format
# ================

module GFF3

import Automa
import Automa.RegExp: @re_str
import BioCore.Exceptions: missingerror
import BioSequences
import BufferedStreams
import GenomicFeatures: GenomicFeatures, Interval, IntervalCollection
import URIParser
importall BioCore

const Bio = BioCore

include("record.jl")
include("reader.jl")
include("writer.jl")

end
