### GFF3

* Reader type: `GenomicFeatures.GFF3.Reader`
* Writer type: `GenomicFeatures.GFF3.Writer`
* Element type: `GenomicFeatures.GFF3.Record`

GFF3 is a text-based file format for representing genomic annotations. The major
difference from BED is that is GFF3 is more structured and can include sequences
in the FASTA file format.

```@docs
GenomicFeatures.GFF3.Reader
GenomicFeatures.GFF3.Writer
GenomicFeatures.GFF3.Record
GenomicFeatures.GFF3.isfeature
GenomicFeatures.GFF3.isdirective
GenomicFeatures.GFF3.iscomment
GenomicFeatures.GFF3.seqid
GenomicFeatures.GFF3.source
GenomicFeatures.GFF3.featuretype
GenomicFeatures.GFF3.seqstart
GenomicFeatures.GFF3.seqend
GenomicFeatures.GFF3.score
GenomicFeatures.GFF3.strand
GenomicFeatures.GFF3.phase
GenomicFeatures.GFF3.attributes
GenomicFeatures.GFF3.content
```
