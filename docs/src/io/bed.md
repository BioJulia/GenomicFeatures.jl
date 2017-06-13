### BED

* Reader type: `BED.Reader`
* Writer type: `BED.Writer`
* Element type: `BED.Record`

BED is a text-based file format for representing genomic annotations like genes,
transcripts, and so on. A BED file has tab-delimited and variable-length fields;
the first three fields denoting a genomic interval are mandatory.

This is an example of RNA transcripts:
```
chr9	68331023	68424451	NM_015110	0	+
chr9	68456943	68486659	NM_001206	0	-
```

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
