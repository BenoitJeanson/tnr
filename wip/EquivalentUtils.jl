using DrWatson
@quickactivate("tnr")

include("../src/GraphUtils.jl")
include("../src/dcpf.jl")
include("Equivalent.jl")
using Logging
using SparseArrays
using MetaGraphsNext
using JuMP
using Gurobi
using Graphs

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

function extract_export_phases(g::MetaGraph, phases::Vector{Float64}, fixed_buses::Vector{String})
    [phases[code_for(g, bus)] for bus in fixed_buses]
end

function convert_id(orig::MetaGraph, target::MetaGraph, vertex::Int)
    code_for(target, label_for(orig, vertex))
end

function convert_id(orig::MetaGraph, target::MetaGraph, vertices::AbstractArray{Int})
    map(v -> convert_id(orig, target, b), vertices)
end

function equivalent_parameters(
    orig::MetaGraph,
    host::MetaGraph,
    buses::Union{AbstractArray{Int},AbstractArray{String}},
    include_branches_limites::Bool;
    outages::Union{AbstractArray{Int},AbstractArray{Tuple{String,String}}}=Int[],)

    _buses_in_orig = isa(buses, AbstractArray{String}) ? to_code(orig, buses) : buses
    _bus_labels = isa(buses, AbstractArray{Int}) ? to_label(orig, buses) : buses
    _buses_in_host = to_code(host, _bus_labels)
    _outages_ids = isa(outages, AbstractArray{Tuple{String,String}}) ?
                   e_code_for(orig, outages) : outages
    _outages_labels = isa(outages, AbstractArray{Int}) ?
                      e_label_for(orig, outages) : outages

    nb_fixed_buses = length(_buses_in_orig)
    flows, _, injections = dc_flow_optim(orig; outages=_outages_ids, fixed_buses=_buses_in_orig, fixed_phases=zeros(Float64, nb_fixed_buses))
    B = zeros(Float64, nb_fixed_buses, nb_fixed_buses)
    F = zeros(Float64, ne(orig), nb_fixed_buses)
    for i in 1:nb_fixed_buses
        dϕ = zeros(Float64, nb_fixed_buses)
        dϕ[i] = 1
        flows2, _, injections2 = dc_flow_optim(orig; outages=_outages_ids, fixed_buses=_buses_in_orig, fixed_phases=dϕ)
        B[:, i] .= injections2 .- injections
        F[:, i] .= flows2 .- flows
    end
    orig_outages = copy(orig)
    foreach(e -> delete!(orig_outages, e[1], e[2]), _outages_labels)
    cc_orig = connected_components(orig_outages)
    ccs = OTS_cc[]
    for cc in cc_orig
        buses = Int[convert_id(orig, host, bus) for bus in _buses_in_orig if bus in cc]
        l = sum([max(orig[label_for(orig, bus)], 0) for bus in cc])
        g = sum([-min(orig[label_for(orig, bus)], 0) for bus in cc])
        push!(ccs, OTS_cc(buses, l, g))
    end
    f_max = [orig[e...].p_max for e in edge_labels(orig)]
    if include_branches_limites
        return Equivalent(_bus_labels, _buses_in_host, injections, B, ccs, flows, f_max, F)
    else
        return Equivalent(_bus_labels, _buses_in_host, injections, B, ccs)
    end
end

