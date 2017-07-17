bigWig
======

Description
-----------

bigWig is a binary file format for associating a floating point number with each
base in the genome. bigWig files are indexed to quickly fetch specific regions.

I/O tools for bigWig are provided from the `GenomicFeatures.BigWig` module,
which exports following three types:
* Reader type: `BigWig.Reader`
* Writer type: `BigWig.Writer`
* Element type: `BigWig.Record`


Examples
--------

A common workflow is to open a file, iterate over records, and close the file:
```julia
# Import the BigWig module.
using GenomicFeatures

# Open a bigWig file (e.g. mapping depth or coverage).
reader = open(BigWig.Reader, "data.cov.bw")

# Iterate over records overlapping with a query interval.
for record in eachoverlap(reader, Interval("Chr2", 5001, 6000))
    # Extract the start position, end position and value of the record,
    startpos = BigWig.chromstart(record)
    endpos = BigWig.chromend(record)
    value = BigWig.value(record)
    # and do something...
end

# Finally, close the reader.
close(reader)
```

`BigWig.values` is a handy function that returns a vector of values. This
returns a value per position within the query region:
```julia
# Get values in Chr2:5001-6000 as a vector of 1000 elements.
BigWig.values(reader, Interval("Chr2", 5001, 6000))
```

Iterating over all records is also supported:
```julia
reader = open(BigWig.Reader, "data.cov.bw")
for record in reader
    # ...
end
close(reader)
```

Creating a bigWig can be written as follows:
```julia
# Open an output file.
file = open("data.cov.bw", "w")

# Initialize a bigWig writer.
writer = BigWig.Writer(file, [("chr1", 2000), ("chr2", 1000)])

# Write records.
write(writer, ("chr1",   1, 100, 1.0))
write(writer, ("chr1", 101, 200, 2.1))
# ...
write(writer, ("chr2",  51, 150, 3.2))

# Close the writer (this closes the file, too).
close(writer)
```


API
---

```@docs
GenomicFeatures.BigWig.Reader
GenomicFeatures.BigWig.chromlist
GenomicFeatures.BigWig.values
GenomicFeatures.BigWig.Writer
GenomicFeatures.BigWig.Record
GenomicFeatures.BigWig.chrom
GenomicFeatures.BigWig.chromid
GenomicFeatures.BigWig.chromstart
GenomicFeatures.BigWig.chromend
GenomicFeatures.BigWig.value
```
