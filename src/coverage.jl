# Coverage
# ========
#
# Coverage of intervals.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

"""
    coverage(intervals)

Compute the coverage of a collection of intervals and return an `IntervalCollection` that contains run-length encoded coverage data.

For example, given intervals like:

    [------]     [------------]
       [---------------]

This function would return a new set of disjoint intervals with annotated coverage like:

    [1][-2-][-1-][--2--][--1--]

# Example

```jldoctest
julia> intervals = [
           Interval("chr1", 1, 8),
           Interval("chr1", 4, 20),
           Interval("chr1", 14, 27)];

julia> coverage(intervals)
IntervalCollection{Interval{UInt32}} with 5 intervals:
  chr1:1-3  .  1
  chr1:4-8  .  2
  chr1:9-13  .  1
  chr1:14-20  .  2
  chr1:21-27  .  1
```
"""
function coverage(stream, seqname_isless::Function=isless)
    cov = IntervalCollection{UInt32}()
    lasts = Int64[]

    stream_next = iterate(stream)
    if stream_next === nothing
        return cov
    end

    current_coverage = 0
    coverage_seqname = ""
    coverage_first = 0
    last_interval_first = 0
    interval, stream_state = stream_next
    stream_next = iterate(stream, stream_state)

    while true
        if seqname(interval) != coverage_seqname
            coverage_process_lasts_heap!(cov, current_coverage, coverage_seqname, coverage_first, lasts)
            if !(isempty(coverage_seqname) || seqname_isless(coverage_seqname, seqname(interval)))
                error("Intervals must be sorted to compute coverage.")
            end

            coverage_seqname = seqname(interval)
            current_coverage = 0
            coverage_first = 0
            last_interval_first = 0
        end

        if leftposition(interval) < last_interval_first
            error("Intervals must be sorted to compute coverage.")
        end

        if !isempty(lasts) && lasts[1] < leftposition(interval)
            pos = DataStructures.heappop!(lasts)
            if leftposition(interval) == pos + 1
                DataStructures.heappush!(lasts, rightposition(interval))
                if stream_next === nothing
                    break
                end
                last_interval_first = leftposition(interval)
                interval, stream_state = stream_next
                stream_next = iterate(stream, stream_state)
            elseif pos == coverage_first - 1
                current_coverage -= 1
            else
                @assert pos >= coverage_first
                push!(cov, Interval{UInt32}(coverage_seqname, coverage_first, pos, STRAND_BOTH, current_coverage))
                current_coverage -= 1
                coverage_first = pos + 1
            end
        else
            if coverage_first == 0
                coverage_first = leftposition(interval)
                current_coverage = 1
            elseif coverage_first == leftposition(interval)
                current_coverage += 1
            else
                if current_coverage > 0
                    push!(cov, Interval{UInt32}(coverage_seqname, coverage_first, leftposition(interval) - 1, STRAND_BOTH, current_coverage))
                end
                current_coverage += 1
                coverage_first = leftposition(interval)
            end

            DataStructures.heappush!(lasts, rightposition(interval))
            if stream_next === nothing
                break
            end
            last_interval_first = leftposition(interval)
            interval, stream_state = stream_next
            stream_next = iterate(stream, stream_state)
        end
    end

    coverage_process_lasts_heap!(cov, current_coverage, coverage_seqname, coverage_first, lasts)

    return cov
end

# Helper function for coverage. Process remaining interval end points after
# all intervals have been read.
function coverage_process_lasts_heap!(cov::IntervalCollection{Interval{UInt32}}, current_coverage, coverage_seqname, coverage_first, lasts)
    while !isempty(lasts)
        pos = DataStructures.heappop!(lasts)
        if pos == coverage_first - 1
            current_coverage -= 1
        else
            @assert pos >= coverage_first
            push!(cov, Interval{UInt32}(coverage_seqname, coverage_first, pos, STRAND_BOTH, current_coverage))
            current_coverage -= 1
            coverage_first = pos + 1
        end
    end
    @assert current_coverage == 0
end
