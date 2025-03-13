using DrWatson
@quickactivate("tnr")

include("../src/GraphUtils.jl")
include("../src/dcpf.jl")
using Logging
using SparseArrays
using MetaGraphsNext
using JuMP
using Gurobi
using ProtoStructs
using Graphs

@proto struct Equivalent
    bus_labels::Vector{String}
    buses::Vector{Int} # ids of the connecting bus when connecting to a target network
    p::Vector{Float64}
    B::Matrix{Float64}
    cc::Vector{Vector{Int}}
end


function subgraph(g::MetaGraph, buses_branches::Union{Dict{Int,Vector{Int}},Dict{String,Vector{Tuple{String,String}}}})

    function _r_add_bus!(h, buses_to_rm, edges_to_rm, g, A, buses_branches, bus)
        bus_label = label_for(g, bus)
        if !haskey(h, bus_label)
            h[bus_label] = g[bus_label]
            !(haskey(buses_branches, bus)) && push!(buses_to_rm, bus_label)
            _edges = haskey(buses_branches, bus) ?
                     buses_branches[bus] :
                     filter(e -> A[bus, e] ≠ 0, 1:size(A)[2])
            for e in _edges
                bus2 = findfirst(i -> i ≠ bus && A[i, e] ≠ 0, 1:size(A)[1])
                _r_add_bus!(h, buses_to_rm, edges_to_rm, g, A, buses_branches, bus2)
                s, d = e_label_for(g, e)
                h[s, d] = g[s, d]
                push!(edges_to_rm, (s, d))
            end
        end
    end

    @info "type of buses_branches", typeof(buses_branches)
    _buses_branches = isa(buses_branches, Dict{Int,Vector{Int}}) ?
                      buses_branches :
                      to_code(g, buses_branches)

    A = incidence_matrix(g; oriented=true)

    sub = _initgraph()
    buses_to_rm = String[]
    edges_to_rm = Tuple{String,String}[]
    _r_add_bus!(sub, buses_to_rm, edges_to_rm, g, A, _buses_branches, first(keys(_buses_branches)))
    _mapedges!(sub)
    remaining = copy(g)
    foreach(e -> delete!(remaining, e[1], e[2]), edges_to_rm)
    foreach(b -> delete!(remaining, b), buses_to_rm)
    _mapedges!(remaining)
    sub, remaining
end

function _graph_to_mat(g::MetaGraph; outages=[], trip=nothing)
    A = _incidence_matrix(g, outages=outages, trip=trip)
    B = [e_index_for(g, e).b for e in edges(g)]
    P = [g[label_for(g, v)] for v in vertices(g)]
    A, B, P
end

function setflows!(g::MetaGraph, flows::Vector{Float64})
    for (i, e) in enumerate(edges(g))
        e_index_for(g, e).p = flows[i]
    end
end

function dc_flow_optim_fixed_phases(A::AbstractMatrix, B::Vector{Float64}, P::Vector{Float64},
    fixed_buses::Vector{Int}, fixed_ϕ::Vector{Float64};
    kwargs...)

    @info "DC PF optim with fixed phases"
    nb_buses, nb_edges = size(A)
    nb_fixed_ϕ = length(fixed_ϕ)
    model = Model(Gurobi.Optimizer)
    # set_silent(model)
    @variable(model, ϕ[1:nb_buses])
    @variable(model, flows[1:nb_edges])
    @constraint(model, [bus in 1:nb_fixed_ϕ], ϕ[fixed_buses[bus]] == fixed_ϕ[bus])
    @constraint(model, [bus in 1:nb_buses; !(bus in fixed_buses)], P[bus] == sum(A[bus, e] * flows[e] for e in 1:nb_edges))
    @constraint(model, [edge in 1:nb_edges], flows[edge] == sum(B[edge] * A[bus, edge] * ϕ[bus] for bus in 1:nb_buses))
    println("MODEL", model)
    optimize!(model)
    if is_solved_and_feasible(model)
        @info "dc_flow_phase solved"
        slack_injections = [sum(A[bus, e] * value(model[:flows][e]) for e in 1:nb_edges) for bus in fixed_buses]
        return value.(model[:flows]), value.(model[:ϕ]), slack_injections
    else
        @warn "dc_flow_phase is not solved"
        return zeros(nb_edges), zeros(nb_buses), zeros(nb_fixed_ϕ)
    end
end

function dc_flow_optim_fixed_phases(g::MetaGraph, fixed_buses::Union{Vector{Int},Vector{String}}, fixed_ϕ::Vector{Float64})
    _fixed_buses = isa(fixed_buses, Vector{String}) ? to_code(g, fixed_buses) : fixed_buses
    A, B, P = _graph_to_mat(g)
    dc_flow_optim_fixed_phases(A, B, P, _fixed_buses, fixed_ϕ)
end

function dc_flow_optim_fixed_phases!(g::MetaGraph, fixed_buses::Union{Vector{Int},Vector{String}}, fixed_ϕ::Vector{Float64})
    flows, ϕ, slack_injections = dc_flow_optim_fixed_phases(g, fixed_buses, fixed_ϕ)
    setflows!(g, flows)
    flows, ϕ, slack_injections
end

function extract_export_phases(g::MetaGraph, phases::Vector{Float64}, fixed_buses::Vector{String})
    [phases[code_for(g, bus)] for bus in fixed_buses]
end

function convert_id(orig::MetaGraph, target::MetaGraph, vertex::Int)
    code_for(target, label_for(orig, vertex))
end

function convert_id(orig::MetaGraph, target::MetaGraph, vertices::AbstractArray{Int})
    map(v -> convert_id(orig, target, b), vertices)
end

function surrogate_parameters(
    orig::MetaGraph,
    host::MetaGraph,
    buses::Union{AbstractArray{Int},AbstractArray{String}};
    outages::Union{AbstractArray{Int},AbstractArray{Tuple{String,String}}}=Int[])

    _buses_in_orig = isa(buses, AbstractArray{String}) ? to_code(orig, buses) : buses
    _bus_labels = isa(buses, AbstractArray{Int}) ? to_label(orig, buses) : buses
    _buses_in_host = to_code(host, _bus_labels)
    _outages_ids = isa(outages, AbstractArray{Tuple{String,String}}) ?
                   e_code_for(orig, outages) : outages
    _outages_labels = isa(outages, AbstractArray{Int}) ?
                      e_label_for(orig, outages) : outages

    nb_fixed_buses = length(_buses_in_orig)
    _, _, injections = dc_flow_optim(orig; outages=_outages_ids, fixed_buses=_buses_in_orig, fixed_phases=zeros(Float64, nb_fixed_buses))
    B = zeros(Float64, nb_fixed_buses, nb_fixed_buses)
    for i in 1:nb_fixed_buses
        dϕ = zeros(Float64, nb_fixed_buses)
        dϕ[i] = 1
        _, _, injections2 = dc_flow_optim(orig; outages=_outages_ids, fixed_buses=_buses_in_orig, fixed_phases=dϕ)
        B[i, :] .= injections2 .- injections
    end
    orig_outages = copy(orig)
    @info "_outages", _outages_ids
    foreach(e -> delete!(orig_outages, e[1], e[2]), _outages_labels)
    cc_orig = connected_components(orig_outages)
    @info "cc_orig", cc_orig
    cc = [Int[convert_id(orig, host, bus) for bus in _buses_in_orig if bus in c] for c in cc_orig]

    Equivalent(_bus_labels, _buses_in_host, injections, B, cc)
end