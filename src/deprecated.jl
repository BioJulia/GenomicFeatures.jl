# Deprecated names
# ================
#
# The squashed names of GenomicFeatures 3.0, which GenomicFeatures names in snake case since 3.1, after the Julia style guide.
# Each is bound to the same function as its snake case name, so methods defined through either name reach callers of both, and every use warns when deprecation warnings are on, as in `Pkg.test`.
# They go in GenomicFeatures 4.0.

# The interval interface, whose functions BioGenerics owns.
Base.@deprecate_binding groupname group_name
Base.@deprecate_binding hasgroupname has_group_name
Base.@deprecate_binding leftposition left_position
Base.@deprecate_binding hasleftposition has_left_position
Base.@deprecate_binding rightposition right_position
Base.@deprecate_binding hasrightposition has_right_position
Base.@deprecate_binding isoverlapping is_overlapping
Base.@deprecate_binding isfilled is_filled

# The legacy names of the group name, which keep BioGenerics' own functions, so that their behaviour doesn't change.
Base.@deprecate_binding seqname BioGenerics.seqname true ", use group_name instead."
Base.@deprecate_binding hasseqname BioGenerics.hasseqname true ", use has_group_name instead."

# The functions that GenomicFeatures owns.
Base.@deprecate_binding eachoverlap each_overlap
Base.@deprecate_binding hasintersection has_intersection
Base.@deprecate_binding intervaltype interval_type false
Base.@deprecate_binding baseintervaltype base_interval_type false
Base.@deprecate_binding metadatatype metadata_type false

# An internal function, which earlier releases' tests and dependents may call qualified.
Base.@deprecate_binding isordered is_ordered false
