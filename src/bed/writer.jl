# BED Writer
# ==========

"""
    BED.Writer(output::IO)

Create a data writer of the BED file format.

# Arguments:
* `output`: data sink
"""
struct Writer <: BioCore.IO.AbstractWriter
    output::IO
end

function BioCore.IO.stream(writer::Writer)
    return writer.output
end

function Base.write(writer::Writer, record::Record)
    return write(writer.output, record, '\n')
end
