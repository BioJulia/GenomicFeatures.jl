# Genomic Interval Manipulation

The `GenomicFeatures` module consists of tools for working efficiently with genomic intervals.

## GenomicInterval Type

Intervals in `GenomicFeatures` are consistent with ranges in Julia: *1-based and end-inclusive*.
When data is read from formats with different representations (i.e. 0-based and/or end-exclusive) they are always converted automatically.
Similarly when writing data, you should not have to reason about off-by-one errors due to format differences while using functionality provided in `GenomicFeatures`.

The [`GenomicInterval`](@ref GenomicInterval) type is defined as
```julia
struct GenomicInterval{T} <: IntervalTrees.AbstractInterval{Int64}
    seqname::String
    first::Int64
    last::Int64
    strand::Strand
    metadata::T
end
```

The first three fields (`seqname`, `first`, and `last`) are mandatory arguments when constructing the [`GenomicInterval`](@ref GenomicInterval) object.
The `seqname` field holds the sequence name associated with the interval.
The `first` and `last` fields are the leftmost and rightmost positions of the interval, which can be accessed with [`leftposition`](@ref leftposition) and [`rightposition`](@ref rightposition) functions, respectively.

The `strand` field can take four kinds of values listed in the next table:

| Symbol | Constant      | Meaning                           |
| :----- | :------------ | :-------------------------------- |
| `'?'`  | `STRAND_NA`   | strandedness is relevant, but unknown |
| `'+'`  | `STRAND_POS`  | positive strand                   |
| `'-'`  | `STRAND_NEG`  | negative strand                   |
| `'.'`  | `STRAND_BOTH` | non-strand-specific feature       |

[`GenomicInterval`](@ref GenomicInterval) is parameterized on metadata type, which lets it efficiently and precisely be specialized to represent intervals from a variety of formats.

The default strand and metadata value are `STRAND_BOTH` and `nothing`:
```jldoctest; setup = :(using GenomicFeatures)
julia> GenomicInterval("chr1", 10000, 20000)
GenomicInterval{Nothing}:
  sequence name: chr1
  leftmost position: 10000
  rightmost position: 20000
  strand: .
  metadata: nothing

julia> GenomicInterval("chr1", 10000, 20000, '+')
GenomicInterval{Nothing}:
  sequence name: chr1
  leftmost position: 10000
  rightmost position: 20000
  strand: +
  metadata: nothing
```

The following example shows all accessor functions for the five fields:
```jldoctest; setup = :(using GenomicFeatures)
julia> i = GenomicInterval("chr1", 10000, 20000, '+', "some annotation")
GenomicInterval{String}:
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


## Collections of GenomicIntervals

Collections of intervals are represented using the [`GenomicIntervalCollection`](@ref GenomicIntervalCollection) type, which is a general purpose indexed container for intervals.
It supports fast intersection operations as well as insertion, deletion, and sorted iteration.

Empty interval collections can be initialized, and intervals elements can be added to the collection one-by-one using `push!`.

```@example
using GenomicFeatures # hide
# The type parameter (Nothing here) indicates the interval metadata type.
col = GenomicIntervalCollection{Nothing}()

for i in 1:100:10000
    push!(col, GenomicInterval("chr1", i, i + 99))
end
```

Incrementally building an interval collection like this works, but [`GenomicIntervalCollection`](@ref GenomicIntervalCollection) also has a bulk insertion constructor that is able to build the indexed data structure extremely efficiently from a sorted vector of intervals.

```jldoctest; setup = :(using GenomicFeatures), output = false
col = GenomicIntervalCollection([Interval("chr1", i, i + 99) for i in 1:100:10000])

# output

GenomicIntervalCollection{GenomicInterval{Nothing}} with 100 intervals:
  chr1:1-100  .  nothing
  chr1:101-200  .  nothing
  chr1:201-300  .  nothing
  chr1:301-400  .  nothing
  chr1:401-500  .  nothing
  chr1:501-600  .  nothing
  chr1:601-700  .  nothing
  chr1:701-800  .  nothing
  ⋮

```

Building [`GenomicIntervalCollection`](@ref GenomicIntervalCollection)s in one shot like this should be preferred when it's convenient or speed is an issue.

## Filtering

Below are some examples of filtering intervals.
The examples take advantage of bulk insertion.
```jldoctest; setup = :(using GenomicFeatures), output = true
intervals = [
  Interval("chr1", 1, 8),
  Interval("chr1", 4, 20),
  Interval("chr1", 14, 27)]

col = IntervalCollection(intervals)

predicate(i) = isodd(leftposition(i))

selected = IntervalCollection(Base.Iterators.filter(predicate, col))
selected = IntervalCollection([x for x in col if predicate(x)])

# output

IntervalCollection{Nothing} with 1 intervals:
  chr1:1-8  .  nothing
```

## Overlap Query

There are number of [`eachoverlap`](@ref eachoverlap) functions in the `GenomicFeatures` module.
They follow two patterns:
- interval versus collection queries which return an iterator over intervals in the collection that overlap the query, and
- collection versus collection queries which iterate over all pairs of overlapping intervals.

```@docs
eachoverlap
```

The order of interval pairs is the same as the following nested loop but [`eachoverlap`](@ref eachoverlap) is often much faster:
```julia
for a in intervals_a, b in intervals_b
    if isoverlapping(a, b)
        # do something...
    end
end
```


## Coverage

A special sort of intersection can also be performed on an interval stream against itself to produce "coverage intervals".

```@docs
coverage
```
