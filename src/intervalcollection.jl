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

mutable struct IntervalCollection{T}
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

    function IntervalCollection{T}() where T
        return new{T}(Dict{String,ICTree{T}}(), 0, ICTree{T}[], false)
    end

    # bulk insertion
    function IntervalCollection{T}(intervals::AbstractVector{Interval{T}}, sort=false) where T
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

function IntervalCollection(intervals::AbstractVector{Interval{T}}, sort=false) where T
    return IntervalCollection{T}(intervals, sort)
end

function IntervalCollection(intervals)
    return IntervalCollection(collect(Interval{metadatatype(intervals)}, intervals), true)
end

function update_ordered_trees!(ic::IntervalCollection{T}) where T
    if ic.ordered_trees_outdated
        ic.ordered_trees = collect(ICTree{T}, values(ic.trees))
        p = sortperm(collect(AbstractString, keys(ic.trees)), lt = isless)
        ic.ordered_trees = ic.ordered_trees[p]
        ic.ordered_trees_outdated = false
    end
end

function Base.push!(ic::IntervalCollection{T}, i::Interval{T}) where T
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

function Base.show(io::IO, ic::IntervalCollection{T}) where T
    n_entries = length(ic)
    println(io, "IntervalCollection{$(T)} with $(n_entries) intervals:")
    if n_entries > 0
        for (k, i) in enumerate(ic)
            if k > 8
                break
            end
            println(IOContext(io, :compact=>true), "  ", i)
        end
        if n_entries > 8
            print(io, "  ⋮")
        end
    end
end

function Base.length(ic::IntervalCollection)
    return ic.length
end

function Base.eltype(::Type{IntervalCollection{T}}) where T
    return Interval{T}
end

function Base.:(==)(a::IntervalCollection{T}, b::IntervalCollection{T}) where T
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

#=
mutable struct IntervalCollectionIteratorState{T}
    i::Int # index into ordered_trees
    tree_state::ICTreeIteratorState{T}

    function IntervalCollectionIteratorState{T}(i::Int) where T
        return new{T}(i)
    end

    function IntervalCollectionIteratorState{T}(i::Int, tree_state) where T
        return new{T}(i, tree_state)
    end
end

function iterinitstate(ic::IntervalCollection{T}) where T
    update_ordered_trees!(ic)
    i = 1
    while i <= length(ic.ordered_trees)
        tree_state = IntervalTrees.iterinitstate(ic.ordered_trees[i])
        if !(tree_state.leaf === nothing || isempty(tree_state.leaf))
            return IntervalCollectionIteratorState{T}(i, tree_state)
        end
        i += 1
    end
    return IntervalCollectionIteratorState{T}(i)
end

function Base.start(ic::IntervalCollection{T}) where T
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
=#

#=
function Base.iterate(ic::IntervalCollection, state=iterinitstate(ic))
    i = state.i
    treeit = iterate(ic.ordered_trees[i], state.tree_state)
    if treeit === nothing
        i += 1
        while i <= length(ic.ordered_trees)
            tree_state = IntervalTrees.iterinitstate(ic.ordered_trees[i])
            if !(tree_state.leaf === nothing || isempty(tree_state.leaf))
                break
            end
            i += 1
        end
    end
    state.i, state.tree_state = i, tree_state
    return value, state
end
=#

function iterprep(ic::IntervalCollection)
    update_ordered_trees!(ic)
    return ()
end

# State is a tuple:
# (ot iteration state, current ot element, current ot element state)
# "ot" is shorthand for "ordered_trees"
# This iterate method basically works like the Iterators.Flatten iterator.
@propagate_inbounds function Base.iterate(ic::IntervalCollection, state = iterprep(ic))
    if state !== ()
        # Iterate over each interval in current ordered tree.
        tree_it = iterate(Base.tail(state)...)
        tree_it !== nothing && return (tree_it[1], (state[1], state[2], tree_it[2]))
    end
    # Iterate over ordered_trees.
    ot_it = (state === () ? iterate(ic.ordered_trees) : iterate(ic.ordered_trees, state[1]))
    ot_it === nothing && return nothing
    iterate(ic, (ot_it[2], ot_it[1]))
end

#=
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
=#

# Filter predicates
# -----------------

true_cmp(a, b) = true

eachoverlap_equal_filter(a, b) = first(a) == first(b) && last(a) == last(b)

# TODO: other prebaked filter predicates


"""
Find a the first interval with matching start and end points.

Returns that interval, or 'nothing' if no interval was found.
"""
function Base.findfirst(a::IntervalCollection{T}, b::Interval{S};
                        filter=true_cmp) where {T,S}
    if !haskey(a.trees, b.seqname)
        return nothing
    else
        return findfirst(a.trees[b.seqname], b, filter)
    end
end


# Overlaps
# --------

function eachoverlap(a::IntervalCollection{T}, b::Interval; filter::F = true_cmp) where {F,T}
    if haskey(a.trees, b.seqname)
        return intersect(a.trees[b.seqname], b)
    else
        return ICTreeIntervalIntersectionIterator{F,T}()
    end
end

function eachoverlap(a::IntervalCollection, b::IntervalCollection; filter = true_cmp)
    seqnames = collect(AbstractString, keys(a.trees) ∩ keys(b.trees))
    sort!(seqnames, lt = isless)
    a_trees = [a.trees[seqname] for seqname in seqnames]
    b_trees = [b.trees[seqname] for seqname in seqnames]
    return IntersectIterator(filter, a_trees, b_trees)
end

struct IntersectIterator{F, S, T}
    filter::F
    a_trees::Vector{ICTree{S}}
    b_trees::Vector{ICTree{T}}
end

#=
mutable struct IntersectIteratorState{F,S,T}
    i::Int  # index into a_trees/b_trees.
    intersect_iterator::ICTreeIntersectionIterator{F,S,T}

    function IntersectIteratorState{F,S,T}(i) where {F,S,T}
        return new{F,S,T}(i)
    end

    function IntersectIteratorState{F,S,T}(i, iter) where {F,S,T}
        return new{F,S,T}(i, iter)
    end
end
=#

function Base.eltype(::Type{IntersectIterator{F,S,T}}) where {F,S,T}
    return Tuple{Interval{S},Interval{T}}
end

function Base.IteratorSize(::Type{IntersectIterator{F,S,T}}) where {F,S,T}
    return Base.SizeUnknown()
end

#=
function Base.start(it::IntersectIterator{F,S,T}) where {F,S,T}
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
=#

# State is a tuple:
# (tree pair iteration state, current intersect iterator, current intersect iterator state)
function Base.iterate(it::IntersectIterator{F, S, T}, state = ()) where {F,S,T}
    if state !== ()
        # Iterate over each intersection between two trees.
        isect = iterate(Base.tail(state)...)
        isect !== nothing && return (isect[1], (state[1], state[2], isect[2]))
    end
    # Iterate over tree pairs.
    treeA = (state === () ? iterate(it.a_trees) : iterate(it.a_trees, state[1]))
    treeB = (state === () ? iterate(it.b_trees) : iterate(it.b_trees, state[1]))
    treeA === nothing && return nothing
    intersect_iterator = intersect(treeA[1], treeB[1], it.filter)
    iterate(it, (treeA[2], intersect_iterator))
end

#=
function Base.next(it::IntersectIterator{F, S, T}, state) where {F,S,T}
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

function Base.done(it::IntersectIterator{F, S, T}, state) where {F, S, T}
    return state.i > length(it.a_trees)
end
=#

function eachoverlap(a, b::IntervalCollection; filter=true_cmp)
    return IntervalCollectionStreamIterator(filter, a, b)
end

struct IntervalCollectionStreamIterator{F,S,T}
    filter::F
    stream::S
    collection::IntervalCollection{T}
end

function Base.eltype(::Type{IntervalCollectionStreamIterator{F,S,T}}) where {F,S,T}
    return Tuple{Interval{metadatatype(S)},Interval{T}}
end

function Base.IteratorSize(::Type{IntervalCollectionStreamIterator{F,S,T}}) where {F,S,T}
    return Base.SizeUnknown()
end

#= ### The old iteration protocol prior to version 0.7 / 1.0 of julia.

mutable struct IntervalCollectionStreamIteratorState{F,Ta,Tb,U}
    intersection::ICTreeIntersection{Tb}
    stream_value::Interval{Ta}
    stream_state::U

    function IntervalCollectionStreamIteratorState{F,Ta,Tb,U}(intersection, stream_value, stream_state) where {F,Ta,Tb,U}
        return new{F,Ta,Tb,U}(intersection, stream_value, stream_state)
    end

    function IntervalCollectionStreamIteratorState{F,Ta,Tb,U}() where {F,Ta,Tb,U}
        return new{F,Ta,Tb,U}(ICTreeIntersection{Tb}())
    end
end

# This mostly follows from SuccessiveTreeIntersectionIterator in IntervalTrees
function Base.start(it::IntervalCollectionStreamIterator{F,S,T}) where {F,S,T}
    stream_state = start(it.stream)
    intersection = ICTreeIntersection{T}()
    while !done(it.stream, stream_state)
        stream_value, stream_state = next(it.stream, stream_state)
        if haskey(it.collection.trees, stream_value.seqname)
            tree = it.collection.trees[stream_value.seqname]
            IntervalTrees.firstintersection!(tree, stream_value, nothing, intersection, it.filter)
            if intersection.index != 0
                return IntervalCollectionStreamIteratorState{F,T,metadatatype(it.stream),typeof(stream_state)}(intersection, stream_value, stream_state)
            end
        end
    end
    return IntervalCollectionStreamIteratorState{S,metadatatype(it.stream),typeof(stream_state)}()
end

function Base.next(it::IntervalCollectionStreamIterator{F,S,T}, state) where {F,S,T}
    intersection = state.intersection
    entry = intersection.node.entries[intersection.index]
    return_value = (state.stream_value, entry)
    IntervalTrees.nextintersection!(intersection.node, intersection.index, state.stream_value, intersection, it.filter)
    while intersection.index == 0 && !done(it.stream, state.stream_state)
        state.stream_value, state.stream_state = next(it.stream, state.stream_state)
        if haskey(it.b.trees, state.stream_value.seqname)
            tree = it.b.trees[state.stream_value.seqname]
            IntervalTrees.firstintersection!(tree, state.stream_value, nothing, intersection, it.filter)
        end
    end
    return return_value, state
end

function Base.done(it::IntervalCollectionStreamIterator, state)
    return state.intersection.index == 0
end
=#

# New julia 0.7 / 1.0 iteration protocol for collection stream iterator.
# State is a tuple: 
# (current_query, stream_state, intersection_object)
function Base.iterate(it::IntervalCollectionStreamIterator{F,S,T}, state = ()) where {F,S,T}
    # If first iteration, make empty intersection, otherwise get it from the state.
    intersection = (state !== () ? state[3] : ICTreeIntersection{T}()) 
    
    # If this is not the first iteration, and there is an available intersection
    # for the current query, return it and search for the next intersection. 
    if state !== () && intersection.index != 0
        entry = intersection.node.entries[intersection.index]
        return_value = (state[1], entry)
        IntervalTrees.nextintersection!(intersection.node, intersection.index, state[1], intersection, it.filter)
        return return_value, (state[1], state[2], intersection)
    end

    # If code reaches this point, there is no valid intersection to return for
    # the current query, so we get the next query and start looking for intersections.
    while intersection.index == 0
        # Get a new query from the stream, and its first intersection in the collection.
        stream_it = (state === () ? iterate(it.stream) : iterate(it.stream, state[2]))
        stream_it === nothing && return nothing # You have reached the end of the stream, stop iterating.
        stream_value = stream_it[1]
        if haskey(it.collection.trees, stream_value.seqname)
            tree = it.collection.trees[stream_value.seqname]
            IntervalTrees.firstintersection!(tree, stream_value, nothing, intersection, it.filter)
            iterate_result = iterate(it, (stream_value, stream_it[2], intersection))
            if iterate_result !== nothing
                return iterate_result
            end
        end
    end
end
