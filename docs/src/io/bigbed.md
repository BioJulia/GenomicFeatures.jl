bigBed
======

Description -----------

bigBed is a binary file format for representing genomic annotations and often
created from BED files. bigBed files are indexed to quickly fetch specific
regions.

I/O tools for bigBed are provided from the `GenomicFeatures.BigBed` module,
which exports following three types:
* Reader type: `BigBed.Reader`
* Writre type: `BigBed.Writer`
* Element type: `BigBed.Record`


Examples
--------

A common workflow is to open a file, iterate over records, and close the file:
```julia
# Import the BigBed module.
using GenomicFeatures

# Open a bigBed file.
reader = open(BigBed.Reader, "data.bb")

# Iterate over records overlapping with a query interval.
for record in eachoverlap(reader, Interval("Chr2", 5001, 6000))
    # Extract the start position, end position and value of the record,
    startpos = BigBed.chromstart(record)
    endpos = BigBed.chromend(record)
    value = BigBed.value(record)
    # and do something...
end

# Finally, close the reader.
close(reader)
```

Iterating over all records is also supported:
```julia
reader = open(BigBed.Reader, "data.bb")
for record in reader
    # ...
end
close(reader)
```

Creating a bigBed file can be done as follows. The `write` call takes a tuple of
3-12 elements (i.e. chromosome name, start position, end position, name, score,
strand, thickstart, thickend, RGB color, blockcount, blocksizes and
blockstarts). The first three are mandatory fields but others are optional.
```julia
# Import RGB type.
using ColorTypes
file = open("data.bb", "w")
writer = BigBed.Writer(file, [("chr1", 1000)])
write(writer, ("chr1", 1, 100, "some name", 100, '+', 10, 90, RGB(0.5, 0.1, 0.2), 2, [4, 10], [10, 20]))
close(writer)
```


API
---

```@docs
GenomicFeatures.BigBed.Reader
GenomicFeatures.BigBed.chromlist
GenomicFeatures.BigBed.Writer
GenomicFeatures.BigBed.Record
GenomicFeatures.BigBed.chrom
GenomicFeatures.BigBed.chromid
GenomicFeatures.BigBed.chromstart
GenomicFeatures.BigBed.chromend
GenomicFeatures.BigBed.name
GenomicFeatures.BigBed.score
GenomicFeatures.BigBed.strand
GenomicFeatures.BigBed.thickstart
GenomicFeatures.BigBed.thickend
GenomicFeatures.BigBed.itemrgb
GenomicFeatures.BigBed.blockcount
GenomicFeatures.BigBed.blocksizes
GenomicFeatures.BigBed.blockstarts
GenomicFeatures.BigBed.optionals
```
