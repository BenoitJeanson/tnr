# module GraphUtils

# export Branch, BranchId, SubBus, BusConf
# export e_index_for, add_subBus, update_edge_ids, build_simple_grid, network2graph, balance!, check_flow_consistency, add_constraint, set_limitations

using Logging
using Graphs
using MetaGraphsNext
using PowerSystems
using Bijections
using ProtoStructs

mutable struct Branch
    b::Float64
    p_max::Float64
    v_nom1::Float64
    v_nom2::Float64
    p::Float64
    outage::Bool
    trip::Bool
end

Branch(b, p_max) = Branch(b, p_max, 0, 0, 0, false, false)
Branch(b, p_max, v_nom1, v_nom2) = Branch(b, p_max, v_nom1, v_nom2, 0, false, false)
Branch(b::Branch) = Branch(b.b, b.p_max, b.v_nom1, b.v_nom2, b.p, false, false)
Branch(b::Branch, p_max) = Branch(b.b, p_max, b.v_nom1, b.v_nom2, b.p, false, false)

struct BranchId
    bus1::Int
    bus2::Int
end

mutable struct SubBus
    p::Float64
    branch_ids::Vector{Int8}
end

mutable struct BusConf
    bus::Int
    subBuses::Vector{SubBus}
end


e_label_for(g::MetaGraph, i::Int) = g[].edges[i]
e_code_for(g::MetaGraph, a::String, b::String) = g[].edges((a, b))

function e_label_for(g::MetaGraph, edge_ids::AbstractArray{Int})
    [e_label_for(g, i) for i in edge_ids]
end

function e_code_for(g::MetaGraph, edges::AbstractArray{Tuple{String,String}})
    [e_code_for(g, e...) for e in edges]
end

function to_label(g::MetaGraph, buses_branches::Dict{Int,<:AbstractArray{Int}})
    dlab = Dict{String,Vector{Tuple{String,String}}}()
    for bus in keys(buses_branches)
        dlab[label_for(g, bus)] = [e_label_for(g, e) for e in buses_branches[bus]]
    end
    dlab
end

function to_code(g::MetaGraph, buses_branches::Dict{String, <:AbstractArray{Tuple{String,String}}})
    dcode = Dict{Int,Vector{Int}}()
    for bus in keys(buses_branches)
        dcode[code_for(g, bus)] = [e_code_for(g, e[1], e[2]) for e in buses_branches[bus]]
    end
    dcode
end

function to_label(g::MetaGraph, buses::AbstractArray{Int})
    [label_for(g, bus) for bus in buses]
end

function to_code(g::MetaGraph, buses::AbstractArray{String})
    [code_for(g, bus) for bus in buses]
end

lsrc(g::MetaGraph, e::Graphs.SimpleEdge) = label_for(g, src(e))
ldst(g::MetaGraph, e::Graphs.SimpleEdge) = label_for(g, dst(e))
e_index_for(g::MetaGraph, e::Graphs.SimpleEdge) = g[lsrc(g, e), ldst(g, e)]

function add_subBus(g::MetaGraph, busconfs::Vector{BusConf})
    h = deepcopy(g)

    e_labels = collect(edge_labels(h))
    edges_origin = collect(edges(h))

    edges_to_change = Dict{Int,Tuple{String,String}}()

    for bc in busconfs
        bus_label = label_for(g, bc.bus)
        for (i, sb) in enumerate(bc.subBuses), br in sb.branch_ids
            new_bus_label = "$bus_label-$i"
            h[bus_label] = g[bus_label] - sb.p
            h[new_bus_label] = sb.p
            elab = get!(edges_to_change, br, e_labels[br])
            edges_to_change[br] = bus_label == elab[1] ? (new_bus_label, elab[2]) : (elab[1], new_bus_label)
            e = edges_origin[br]
            rem_edge!(h, src(e), dst(e))
        end
    end

    for (e, elab) in edges_to_change
        h[elab...] = g[e_labels[e]...]
    end

    h
end

# TODO: consider using e_label_for 
function update_edge_ids(edge_ids::Vector{Int}, id::Int, g_orig::MetaGraph, g_modified::MetaGraph)
    function extract_to_dash(st::String)
        fl = findlast('-', st)
        isnothing(fl) ? st : st[1:fl-1]
    end

    if id ≠ 0
        push!(edge_ids, id)
    end
    result = zeros(Int, length(edge_ids))

    edges_o = collect(edge_labels(g_orig))
    edges_m = [(extract_to_dash(el[1]), extract_to_dash(el[2])) for el in edge_labels(g_modified)]

    ids_labels = [edges_o[id] for id in edge_ids]

    for (i_m, e_modif) in enumerate(edges_m)
        for (i_o, elo) in enumerate(ids_labels)
            if elo == e_modif
                result[i_o] = i_m
                continue
            end
        end
    end
    if id ≠ 0
        id_res = pop!(result)
    else
        id_res = 0
    end
    result, id_res
end

function update_edge_ids(edge_ids::Vector{Int}, g_orig::MetaGraph, g_modified::MetaGraph)
    update_edge_ids(edge_ids, 0, g_orig, g_modified)[1]
end


function add_subBus(g::MetaGraph, busconfs::Vector{BusConf}, edge_ids::Vector{Int})
    h = add_subBus(g, busconfs)
    h, update_edge_ids(edge_ids, g, h)
end

function add_subBus(g::MetaGraph, busconfs::Vector{BusConf}, edge_ids::Vector{Int}, id::Int)
    h = add_subBus(g, busconfs)
    uei = update_edge_ids(edge_ids, id, g, h)
    h, uei[1], uei[2]
end

function _initgraph()
    MetaGraph(
        DiGraph();
        label_type=String,
        vertex_data_type=Float64,
        edge_data_type=Branch,
        graph_data=(
            edges=Bijection{Int,Tuple{String,String}}(),)
    )
end

function _mapedges!(g::MetaGraph)
    gd = g[].edges
    foreach(k -> delete!(gd, k), keys(gd))
    foreach(ie -> gd[ie[1]] = ie[2], enumerate(edge_labels(g)))
end


function build_simple_grid(; micro=true)
    g = _initgraph()
    g["1"] = micro ? -2 : -3
    for i = 2:(micro ? 3 : 4)
        g["$i"] = 1
    end
    g["1", "2"] = Branch(1, 1)
    g["2", "3"] = Branch(1, 1)
    if micro
        g["1", "3"] = Branch(1, 1)
    else
        g["3", "4"] = Branch(1, 1)
        g["2", "4"] = Branch(1, 1)
        g["1", "4"] = Branch(1, 1)
    end
    _mapedges!(g)
    g
end

function network2graph(sys::System, buslabels=num -> "$num")
    g = _initgraph()

    init_p_bus(g, sys, objType, op) =
        for k in get_components(objType, sys)
            label = buslabels(k.bus.number)
            @debug "bus $label\t$(repr(op)) $(k.active_power)"
            g[label] = op(haskey(g, label) ? g[label] : 0, k.active_power)
        end
    init_p_bus(g, sys, PowerLoad, +)
    init_p_bus(g, sys, Generator, -)

    @debug "branches"
    for k in get_components(ACBranch, sys)
        labelfrom = buslabels(k.arc.from.number)
        labelto = buslabels(k.arc.to.number)
        if !haskey(g, labelfrom)
            @debug "bus: $label"
            g[labelfrom] = 0
        end
        if !haskey(g, labelto)
            @debug "bus: $labelto"
            g[labelto] = 0
        end
        @debug "(\"$labelfrom\", \"$labelto\") "
        if labelfrom < labelto
            g[labelfrom, labelto] = Branch(1 / k.x, 10, k.arc.from.base_voltage, k.arc.to.base_voltage)
        else
            g[labelto, labelfrom] = Branch(1 / k.x, 10, k.arc.to.base_voltage, k.arc.from.base_voltage)
        end
    end

    to_remove = []
    for node in vertices(g)
        if isempty(inneighbors(g, node)) && isempty(outneighbors(g, node))
            push!(to_remove, label_for(g, node))
        end
    end
    @debug !isempty(to_remove) && "\nremove nodes: $to_remove"

    foreach(node -> rem_vertex!(g, code_for(g, node)), to_remove)

    _mapedges!(g)
    g
end

function edge_label_to_num(g, labsrc, labdst)
    d = getindex(g)
    haskey(d, (labsrc, labdst)) ? d[(labsrc, labdst)] : nothing
end

@enum BalanceType all_non_zero_uniform gen_proportional

function balance!(g, btype::BalanceType=gen_proportional)
    non_zeros = [label for label in labels(g) if g[label] ≠ 0]
    imbalance = sum(g[label] for label in non_zeros)

    if btype == all_non_zero_uniform
        uniform = imbalance / length(non_zeros)
        for label in non_zeros
            g[label] -= uniform
        end

    elseif btype == gen_proportional
        generators = [label for label in labels(g) if g[label] ≤ 0]
        if isempty(generators)
            println("no generator to balance")
        else
            total_gen = sum(g[label] for label in generators)
            for label in generators
                g[label] -= imbalance * g[label] / total_gen
            end
        end

    else
        println("BALANCE TYPE TO BE IMPLEMENTED")
    end #TODO oter types
end


function check_flow_consistency(g; v::Bool=false)
    flows = Dict(l => g[l...].p for l in edge_labels(g))
    injections = Dict(l => g[l] for l in labels(g))

    # Initialize net flows for each node
    net_flows = Dict{String,Float64}()

    # Update net flows based on branch flows
    for ((from, to), flow) in flows
        net_flows[from] = get(net_flows, from, 0.0) - flow
        net_flows[to] = get(net_flows, to, 0.0) + flow
    end

    # Compare with injections
    consistent = true

    v && println("\nComparison with injections:")
    for (node, injection) in injections
        net_flow = get(net_flows, node, 0.0)

        v && println("Node $node: Injection = $injection, Net Flow = $net_flow, Difference = $(injection - net_flow)")
        if abs(injection - net_flow) > 1e-9  # Tolerance for floating-point comparison
            consistent = false
        end
    end

    if v
        if consistent
            println("\nThe flows and injections are consistent.")
        else
            println("\nThe flows and injections are not consistent.")
        end
    end
    consistent
end

add_constraint(g, fun) = foreach(e -> fun(e_index_for(g, e)), edges(g))

set_limitations(branches, lim) = Dict(k => begin
    br_out = Branch(br_in, br_in.v_nom1 == 138.0 ? lim : 10)
end for (k, br_in) in branches)


function connected(g::MetaGraph, bus_origin::String; outages=Int[], trip=nothing)

    """
    _extractsubgraph is embedded as it only applies if the bus_labels are the labels of buses belonging to a single connected component 
    """
    function _extractsubgraph(g::MetaGraph, bus_labels::Vector{String})
        g_sub = _initgraph()
        for bus in bus_labels
            g_sub[bus] = g[bus]
        end
        for e in edges(g)
            slabel = label_for(g, src(e))
            if slabel in bus_labels
                dlabel = label_for(g, dst(e))
                g_sub[slabel, dlabel] = g[slabel, dlabel]
            end
        end
        g_sub
    end

    g_result = copy(g)

    openbranches = copy(outages)
    if !isnothing(trip)
        push!(openbranches, trip)
    end

    for (i, e) in [(i, e) for (i, e) in enumerate(edges(g_result)) if i in openbranches]
        rem_edge!(g_result, src(e), dst(e))
    end

    ccs = connected_components(g_result)
    cc = ccs[findfirst(cc -> bus_origin in [label_for(g_result, v) for v in cc], ccs)]
    _extractsubgraph(g_result, [label_for(g_result, v) for v in cc])
end

function total_load(g)
    result = 0
    for v in vertices(g)
        lab = label_for(g, v)
        if g[lab] > 0
            result += g[lab]
        end
    end
    result
end
