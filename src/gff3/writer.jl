# GFF3 Writer
# ===========

"""
    GFF3.Writer(output::IO)

Create a data writer of the GFF3 file format.

# Arguments:
* `output`: data sink
"""
mutable struct Writer{T<:IO} <: BioCore.IO.AbstractReader
    stream::T
end

function BioCore.IO.stream(writer::Writer)
    return writer.stream
end

function Base.write(writer::Writer, record::Record)
    return write(writer.stream, record, '\n')
end
