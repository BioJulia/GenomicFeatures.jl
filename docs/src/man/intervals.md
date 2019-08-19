# Genomic Interval Manipulation


The `GenomicFeatures` module consists of tools for working efficiently with
genomic intervals.


## Interval Type

Intervals in GenomicFeatures.jl are consistent with ranges in Julia: *1-based
and end-inclusive*. When data is read from formats with different
representations (i.e. 0-based and/or end-exclusive) they are always converted
automatically.  Similarly when writing data. You should not have to reason about
off-by-one errors due to format differences while using functionality provided
in GenomicFeatures.jl.

The `Interval` type is defined as
```julia
struct Interval{T} <: AbstractInterval{Int64}
    seqname::String
    first::Int64
    last::Int64
    strand::Strand
    metadata::T
end
```

The first three fields (`seqname`, `first`, and `last`) are mandatory arguments
when constructing an `Interval` object. `seqname` is the sequence name
associated with the interval. The `first` and `last` fields are the leftmost and
rightmost positions of the interval, which can be accessed with `leftposition`
and `rightposition` functions, respectively.

The `strand` field can take four kinds of values listed in the next table:

| Symbol | Constant      | Meaning                           |
| :----- | :------------ | :-------------------------------- |
| `'?'`  | `STRAND_NA`   | strand is unknown or inapplicable |
| `'+'`  | `STRAND_POS`  | positive strand                   |
| `'-'`  | `STRAND_NEG`  | negative strand                   |
| `'.'`  | `STRAND_BOTH` | non-strand-specific feature       |

`Interval` is parameterized on metadata type, which lets it efficiently and
precisely be specialized to represent intervals from a variety of formats.

The default strand and metadata value are `STRAND_BOTH` and `nothing`:
```jlcon
julia> Interval("chr1", 10000, 20000)
GenomicFeatures.Interval{Void}:
  sequence name: chr1
  leftmost position: 10000
  rightmost position: 20000
  strand: .
  metadata: nothing

julia> Interval("chr1", 10000, 20000, '+')
GenomicFeatures.Interval{Void}:
  sequence name: chr1
  leftmost position: 10000
  rightmost position: 20000
  strand: +
  metadata: nothing

```

The following example shows all accessor functions for the five fields:
```jlcon
julia> i = Interval("chr1", 10000, 20000, '+', "some annotation")
GenomicFeatures.Interval{String}:
  sequence name: chr1
  leftmost position: 10000
  rightmost position: 20000
  strand: +
  metadata: some annotation

julia> seqname(i)
"chr1"

julia> leftposition(i)
10000

julia> rightposition(i)
20000

julia> strand(i)
STRAND_POS

julia> metadata(i)
"some annotation"

```


## Collections of Intervals

Collections of intervals are represented using the `IntervalCollection` type,
which is a general purpose indexed container for intervals. It supports fast
intersection operations as well as insertion, deletion, and sorted iteration.

Interval collections can be initialized by inserting elements one by one using
`push!`.

```julia
# The type parameter (Void here) indicates the interval metadata type.
col = IntervalCollection{Void}()

for i in 1:100:10000
    push!(col, Interval("chr1", i, i + 99))
end
```

Incrementally building an interval collection like this works, but
`IntervalCollection` also has a bulk insertion constructor that is able to build
the indexed data structure extremely efficiently from an array of intervals.

```julia
col = IntervalCollection([Interval("chr1", i, i + 99) for i in 1:100:10000])
```

Bulding `IntervalCollections` in one shot like this should be preferred when
it's convenient or speed in an issue.

`IntervalCollection`s can also be build directly from a genome annotation file,
here in GFF3 format:

```julia
reader = open(GFF3.Reader, "some_genome.gff3")
features = IntervalCollection(reader)
```


## Overlap Query

There are number of `eachoverlap` function in the `GenomicFeatures` module. They
follow two patterns: interval versus collection queries which return an iterator
over intervals in the collection that overlap the query, and collection versus
collection queries which iterate over all pairs of overlapping intervals.

```@docs
eachoverlap
```

The order of interval pairs is the same as the following nested loop but
`eachoverlap` is often much faster:
```julia
for a in intervals_a, b in intervals_b
    if isoverlapping(a, b)
        # do something...
    end
end
```


## Coverage

A special sort of intersection can also be performed on an interval stream
against itself to produce "coverage intervals".

```@docs
coverage
```
