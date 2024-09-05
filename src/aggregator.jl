"""
    abstract type Aggregator

Abstract type for aggregators of PTMC simulation results.
Concrete types must dispatch the `create_aggregation` and `update_aggregation` functions.
See `HistoryAggregator` for an example.
"""
abstract type Aggregator end

"""
    struct HistoryAggregator <: Aggregator

An aggregator that stores the history of each photon in the simulation.
"""
struct HistoryAggregator <: Aggregator
    history::Vector{Vector{Photon}}
    HistoryAggregator() = new([])
end

create_aggregation(agg::HistoryAggregator, pstart::Photon) = push!(agg.history, [pstart])
update_aggregation(agg::HistoryAggregator, pstart::Photon, pend::Photon) = push!(agg.history[end], pend)