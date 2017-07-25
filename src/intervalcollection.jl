# An IntervalCollection is an efficiently stored and indexed set of annotated
# genomic intervals. It looks something like this.
#
#                                      ┌─────┐
#                                      │trees│
#                                      └─────┘
#                                         │
#                              ┌──────────┴──────────┬────────────┐
#                              ▼                     ▼            │
#  Each sequence has       ┌──────┐              ┌──────┐         ▼
#    an associated         │ chr1 │              │ chr2 │
#    IntervalTree     ┌────┴──────┴────┐    ┌────┴──────┴────┐    ...
#                     │ IntervalTree 1 │    │ IntervalTree 2 │
#                     └────────────────┘    └────────────────┘
#                              │
#                    ┌─────────┴─────────────────────────┬────────────────────────┐
#                    │                                   │                        │
#                    ▼                                   ▼                        ▼
#   ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓  ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#   ┃ Interval{T}(10000, 20000, ...) ┃  ┃ Interval{T}(35000, 40000, ...) ┃      ...
#   ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛  ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
#
#
#
#
#                               ┌────────────────┐
#       ordered_trees holds     │ IntervalTree 1 │
#         an array of the       ├────────────────┤
#       some IntervalTrees,     │ IntervalTree 2 │
#           ordered by          └────────────────┘
#       chromosome for fast
#       ordered iteration.             ...
#

# Aliases for types of IntervalTrees.jl (IC: Interval Collection).
const ICTree{T}                               = IntervalTrees.IntervalBTree{Int64,Interval{T},64}
const ICTreeIteratorState{T}                  = IntervalTrees.IntervalBTreeIteratorState{Int64,Interval{T},64}
const ICTreeIntersection{T}                   = IntervalTrees.Intersection{Int64,Interval{T},64}
const ICTreeIntersectionIterator{F,S,T}       = IntervalTrees.IntersectionIterator{F,Int64,Interval{S},64,Interval{T},64}
const ICTreeIntervalIntersectionIterator{F,T} = IntervalTrees.IntervalIntersectionIterator{F, Int64,Interval{T},64}

type IntervalCollection{T}
    # Sequence name mapped to IntervalTree, which in turn maps intervals to
    # a list of metadata.
    trees::Dict{String,ICTree{T}}

    # Keep track of the number of stored intervals
    length::Int

    # A vector of values(trees) sorted on sequence name.
    # This is used to iterate intervals as efficiently as possible, but is only
    # updated as needed, indicated by the ordered_trees_outdated flag.
    ordered_trees::Vector{ICTree{T}}
    ordered_trees_outdated::Bool

    function (::Type{IntervalCollection{T}}){T}()
        return new{T}(Dict{String,ICTree{T}}(), 0, ICTree{T}[], false)
    end

    # bulk insertion
    function (::Type{IntervalCollection{T}}){T}(intervals::AbstractVector{Interval{T}}, sort=false)
        if sort
            sort!(intervals)
        else
            if !issorted(intervals)
                error("Intervals must be sorted, or `sort=true` set, to construct an IntervalCollection")
            end
        end

        n = length(intervals)
        trees = Dict{String,ICTree{T}}()
        i = 1
        while i <= n
            j = i
            while j <= n && intervals[i].seqname == intervals[j].seqname
                j += 1
            end
            trees[intervals[i].seqname] = ICTree{T}(view(intervals, i:j-1))
            i = j
        end
        return new{T}(trees, n, ICTree{T}[], true)
    end
end

function IntervalCollection{T}(intervals::AbstractVector{Interval{T}}, sort=false)
    return IntervalCollection{T}(intervals, sort)
end

function IntervalCollection(intervals)
    return IntervalCollection(collect(Interval{metadatatype(intervals)}, intervals), true)
end

function update_ordered_trees!{T}(ic::IntervalCollection{T})
    if ic.ordered_trees_outdated
        ic.ordered_trees = collect(ICTree{T}, values(ic.trees))
        p = sortperm(collect(AbstractString, keys(ic.trees)), lt=isless)
        ic.ordered_trees = ic.ordered_trees[p]
        ic.ordered_trees_outdated = false
    end
end

function Base.push!{T}(ic::IntervalCollection{T}, i::Interval{T})
    if !haskey(ic.trees, i.seqname)
        tree = ICTree{T}()
        ic.trees[i.seqname] = tree
        ic.ordered_trees_outdated = true
    else
        tree = ic.trees[i.seqname]
    end
    push!(tree, i)
    ic.length += 1
    return ic
end

function Base.show{T}(io::IO, ic::IntervalCollection{T})
    n_entries = length(ic)
    println(io, "IntervalCollection{$(T)} with $(n_entries) intervals:")
    if n_entries > 0
        for (k, i) in enumerate(ic)
            if k > 8
                break
            end
            println(IOContext(io, compact=true), "  ", i)
        end
        if n_entries > 8
            print(io, "  ⋮")
        end
    end
end

function Base.length(ic::IntervalCollection)
    return ic.length
end

function Base.eltype{T}(::Type{IntervalCollection{T}})
    return Interval{T}
end

function Base.:(==){T}(a::IntervalCollection{T}, b::IntervalCollection{T})
    if length(a) != length(b)
        return false
    end
    for (i, j) in zip(a, b)
        if i != j
            return false
        end
    end
    return true
end


# Iterators
# ---------

type IntervalCollectionIteratorState{T}
    i::Int # index into ordered_trees
    tree_state::ICTreeIteratorState{T}

    function (::Type{IntervalCollectionIteratorState{T}}){T}(i::Int)
        return new{T}(i)
    end

    function (::Type{IntervalCollectionIteratorState{T}}){T}(i::Int, tree_state)
        return new{T}(i, tree_state)
    end
end

function Base.start{T}(ic::IntervalCollection{T})
    update_ordered_trees!(ic)
    i = 1
    while i <= length(ic.ordered_trees)
        tree_state = start(ic.ordered_trees[i])
        if !done(ic.ordered_trees[i], tree_state)
            return IntervalCollectionIteratorState{T}(i, tree_state)
        end
        i += 1
    end
    return IntervalCollectionIteratorState{T}(i)
end

function Base.next(ic::IntervalCollection, state)
    i = state.i
    value, tree_state = next(ic.ordered_trees[i], state.tree_state)
    if done(ic.ordered_trees[i], tree_state)
        i += 1
        while i <= length(ic.ordered_trees)
            tree_state = start(ic.ordered_trees[i])
            if !done(ic.ordered_trees[i], tree_state)
                break
            end
            i += 1
        end
    end
    state.i, state.tree_state = i, tree_state
    return value, state
end

function Base.done(ic::IntervalCollection, state)
    return state.i > length(ic.ordered_trees)
end


# Filter predicates
# -----------------

true_cmp(a, b) = true

eachoverlap_equal_filter(a, b) = first(a) == first(b) && last(a) == last(b)

# TODO: other prebaked filter predicates


"""
Find a the first interval with matching start and end points.

Returns that interval wrapped in a Nullable, or an empty Nullable if no such
interval was found.
"""
function Base.findfirst{T,S}(a::IntervalCollection{T}, b::Interval{S};
                             filter=true_cmp)
    if !haskey(a.trees, b.seqname)
        return Nullable{Interval{T}}()
    else
        return findfirst(a.trees[b.seqname], b, filter)
    end
end


# Overlaps
# --------

function eachoverlap{F,T}(a::IntervalCollection{T}, b::Interval; filter::F=true_cmp)
    if haskey(a.trees, b.seqname)
        return intersect(a.trees[b.seqname], b)
    else
        return ICTreeIntervalIntersectionIterator{F,T}()
    end
end

function eachoverlap(a::IntervalCollection, b::IntervalCollection; filter=true_cmp)
    seqnames = collect(AbstractString, keys(a.trees) ∩ keys(b.trees))
    sort!(seqnames, lt=isless)
    a_trees = [a.trees[seqname] for seqname in seqnames]
    b_trees = [b.trees[seqname] for seqname in seqnames]
    return IntersectIterator(filter, a_trees, b_trees)
end

immutable IntersectIterator{F, S, T}
    filter::F
    a_trees::Vector{ICTree{S}}
    b_trees::Vector{ICTree{T}}
end

type IntersectIteratorState{F,S,T}
    i::Int  # index into a_trees/b_trees.
    intersect_iterator::ICTreeIntersectionIterator{F,S,T}

    function (::Type{IntersectIteratorState{F,S,T}}){F,S,T}(i)
        return new{F,S,T}(i)
    end

    function (::Type{IntersectIteratorState{F,S,T}}){F,S,T}(i, iter)
        return new{F,S,T}(i, iter)
    end
end

function Base.eltype{F,S,T}(::Type{IntersectIterator{F,S,T}})
    return Tuple{Interval{S},Interval{T}}
end

function Base.iteratorsize{F,S,T}(::Type{IntersectIterator{F,S,T}})
    return Base.SizeUnknown()
end

function Base.start{F,S,T}(it::IntersectIterator{F,S,T})
    i = 1
    while i <= length(it.a_trees)
        intersect_iterator = intersect(it.a_trees[i], it.b_trees[i], it.filter)
        intersect_iterator_state = start(intersect_iterator)
        if !done(intersect_iterator, intersect_iterator_state)
            return IntersectIteratorState{F,S,T}(i, intersect_iterator)
        end
        i += 1
    end
    return IntersectIteratorState{F,S,T}(i)
end

function Base.next{F,S,T}(it::IntersectIterator{F, S, T}, state)
    i, intersect_iterator = state.i, state.intersect_iterator
    value, intersect_iterator_state = next(intersect_iterator, nothing)
    if done(intersect_iterator, intersect_iterator_state)
        i += 1
        while i <= length(it.a_trees)
            intersect_iterator = intersect(it.a_trees[i], it.b_trees[i],
                                           it.filter)
            intersect_iterator_state = start(intersect_iterator)
            if !done(intersect_iterator, intersect_iterator_state)
                break
            end
            i += 1
        end
    end
    state.i, state.intersect_iterator = i, intersect_iterator
    return value, state
end

function Base.done{F, S, T}(it::IntersectIterator{F, S, T}, state)
    return state.i > length(it.a_trees)
end

function eachoverlap(a, b::IntervalCollection; filter=true_cmp)
    return IntervalCollectionStreamIterator(filter, a, b)
end

immutable IntervalCollectionStreamIterator{F,S,T}
    filter::F
    a::S
    b::IntervalCollection{T}
end

function Base.eltype{F,S,T}(::Type{IntervalCollectionStreamIterator{F,S,T}})
    return Tuple{Interval{metadatatype(S)},Interval{T}}
end

function Base.iteratorsize{F,S,T}(::Type{IntervalCollectionStreamIterator{F,S,T}})
    return Base.SizeUnknown()
end

type IntervalCollectionStreamIteratorState{F,Ta,Tb,U}
    intersection::ICTreeIntersection{Tb}
    a_value::Interval{Ta}
    a_state::U

    function (::Type{IntervalCollectionStreamIteratorState{F,Ta,Tb,U}}){F,Ta,Tb,U}(intersection, a_value, a_state)
        return new{F,Ta,Tb,U}(intersection, a_value, a_state)
    end

    function (::Type{IntervalCollectionStreamIteratorState{F,Ta,Tb,U}}){F,Ta,Tb,U}()
        return new{F,Ta,Tb,U}(ICTreeIntersection{Tb}())
    end
end

# This mostly follows from SuccessiveTreeIntersectionIterator in IntervalTrees
function Base.start{F,S,T}(it::IntervalCollectionStreamIterator{F,S,T})
    a_state = start(it.a)
    intersection = ICTreeIntersection{T}()
    while !done(it.a, a_state)
        a_value, a_state = next(it.a, a_state)
        if haskey(it.b.trees, a_value.seqname)
            tree = it.b.trees[a_value.seqname]
            IntervalTrees.firstintersection!(tree, a_value, Nullable{Interval{T}}(), intersection, it.filter)
            if intersection.index != 0
                return IntervalCollectionStreamIteratorState{F,T,metadatatype(it.a),typeof(a_state)}(intersection, a_value, a_state)
            end
        end
    end
    return IntervalCollectionStreamIteratorState{S,metadatatype(it.a),typeof(a_state)}()
end

function Base.next{F,S,T}(it::IntervalCollectionStreamIterator{F,S,T}, state)
    intersection = state.intersection
    entry = intersection.node.entries[intersection.index]
    return_value = (state.a_value, entry)
    IntervalTrees.nextintersection!(intersection.node, intersection.index, state.a_value, intersection, it.filter)
    while intersection.index == 0 && !done(it.a, state.a_state)
        state.a_value, state.a_state = next(it.a, state.a_state)
        if haskey(it.b.trees, state.a_value.seqname)
            tree = it.b.trees[state.a_value.seqname]
            IntervalTrees.firstintersection!(tree, state.a_value, Nullable{Interval{T}}(), intersection, it.filter)
        end
    end
    return return_value, state
end

function Base.done(it::IntervalCollectionStreamIterator, state)
    return state.intersection.index == 0
end
