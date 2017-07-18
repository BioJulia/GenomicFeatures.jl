BED
===

Description
-----------

BED is a text-based file format for representing genomic annotations like genes,
transcripts, and so on. A BED file has tab-delimited and variable-length fields;
the first three fields denoting a genomic interval are mandatory.

This is an example of RNA transcripts:
```
chr9	68331023	68424451	NM_015110	0	+
chr9	68456943	68486659	NM_001206	0	-
```

I/O tools for BED are provided from the `GenomicFeatures.BED` module,
which exports following three types:
* Reader type: `BED.Reader`
* Writer type: `BED.Writer`
* Element type: `BED.Record`


Examples
--------

Here is a common workflow to iterate over all records in a BED file:
```julia
# Import the BED module.
using GenomicFeatures

# Open a BED file.
reader = open(BED.Reader, "data.bed")

# Iterate over records.
for record in reader
    # Do something on record (see Accessors section).
    chrom = BED.chrom(record)
    # ...
end

# Finally, close the reader.
close(reader)
```

If you repeatedly access records within specific ranges, it would be more
efficient to construct an `IntervalCollection` object from a BED reader:
```julia
# Create an interval collection in memory.
icol = open(BED.Reader, "data.bed") do reader
    IntervalCollection(reader)
end

# Query overlapping records.
for interval in eachoverlap(icol, Interval("chrX", 40001, 51500))
    # A record is stored in the metadata field of an interval.
    record = metadata(interval)
    # ...
end
```


API
---

```@docs
GenomicFeatures.BED.Reader
GenomicFeatures.BED.Writer
GenomicFeatures.BED.Record
GenomicFeatures.BED.chrom
GenomicFeatures.BED.chromstart
GenomicFeatures.BED.chromend
GenomicFeatures.BED.name
GenomicFeatures.BED.score
GenomicFeatures.BED.strand
GenomicFeatures.BED.thickstart
GenomicFeatures.BED.thickend
GenomicFeatures.BED.itemrgb
GenomicFeatures.BED.blockcount
GenomicFeatures.BED.blocksizes
GenomicFeatures.BED.blockstarts
```
