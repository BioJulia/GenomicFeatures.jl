# Snake case names
# ================
#
# Julia's style guide asks packages to write function names in snake case.
# The interval interface belongs to BioGenerics, whose functions GenomicFeatures extends; 3.0 exported them under their squashed names, which src/deprecated.jl deprecates.
# Each snake case name below is bound to the BioGenerics function itself, so the two names are one function: a method defined through either name reaches callers of both.

"""
    group_name(x)

Get the name of the group of `x`, such as the chromosome of an interval.
It is BioGenerics' `groupname`, so a method defined through either name reaches callers of both, though GenomicFeatures deprecates the name `groupname`.
"""
const group_name = BioGenerics.groupname

"""
    has_group_name(x)::Bool

Test if `x` has a group name.
It is BioGenerics' `hasgroupname`, whose name GenomicFeatures deprecates.
"""
const has_group_name = BioGenerics.hasgroupname

"""
    left_position(x)

Get the leftmost position of `x`, 1-based.
It is BioGenerics' `leftposition`, whose name GenomicFeatures deprecates.
"""
const left_position = BioGenerics.leftposition

"""
    has_left_position(x)::Bool

Test if `x` has a leftmost position.
It is BioGenerics' `hasleftposition`, whose name GenomicFeatures deprecates.
"""
const has_left_position = BioGenerics.hasleftposition

"""
    right_position(x)

Get the rightmost position of `x`, 1-based and inclusive.
It is BioGenerics' `rightposition`, whose name GenomicFeatures deprecates.
"""
const right_position = BioGenerics.rightposition

"""
    has_right_position(x)::Bool

Test if `x` has a rightmost position.
It is BioGenerics' `hasrightposition`, whose name GenomicFeatures deprecates.
"""
const has_right_position = BioGenerics.hasrightposition

"""
    is_overlapping(a, b)::Bool

Test if `a` and `b` overlap.
It is BioGenerics' `isoverlapping`, whose name GenomicFeatures deprecates.
"""
const is_overlapping = BioGenerics.isoverlapping

"""
    is_filled(x)::Bool

Test if `x` holds data, rather than being empty.
It is BioGenerics' `isfilled`, whose name GenomicFeatures deprecates.
"""
const is_filled = BioGenerics.isfilled
