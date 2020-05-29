# Filter
# ======
#
# Filter function for IntervalCollection
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

"""
    filter(f, coll:IntervalCollection)

Return a copy of `coll`, removing elements (Intervals) for which `f` is false. 
The function `f` is passed one argument (Interval).

This is an eager implimentation of filtering for IntervalCollection, 
which is equivalent to the following code:

```julia
result=IntervalCollection{T}()   # T is the type of metadata of the element in coll
for i in Base.Iterators.filter(f, coll)
    push!(result,i)
end
result
```

# Examples
# --------

This function can be used with `coverage` function
to select intervals with specified number of coverage.

The following code selects non-overlapping intervals.

```julia
coll=IntervalCollection{Nothing}()
push!(coll, Interval("chr1", 1, 9))
push!(coll, Interval("chr1", 4, 7))
cov=coverage(coll)
cov1=filter(cov) do i
    i.metadata==1  # specify coverage
end
@show cov1
# =>
# IntervalCollection{UInt32} with 2 intervals:
#  chr1:1-3  .  1
#  chr1:8-9  .  1
```
"""

function Base.filter(f::Function, coll::IntervalCollection{T}) where T
    result=IntervalCollection{T}()
    for c in coll
        if f(c)
            push!(result,c)
        end
    end
    result
end

