include("../src/GraphUtils.jl")
include("../src/dcpf.jl")
include("Equivalent.jl")
using Logging
using SparseArrays
using MetaGraphsNext
using JuMP
using Gurobi
using Graphs

function subgraph(g::MetaGraph, buses_branches::Dict{VLabel,Vector{ELabel}})

    function _r_add_bus!(h, buses_to_rm, edges_to_rm, g, _buses_branches, bus_label)
        if !haskey(h, bus_label)
            h[bus_label] = g[bus_label]
            !(haskey(_buses_branches, bus_label)) && push!(buses_to_rm, bus_label)
            _edges = haskey(_buses_branches, bus_label) ?
                     _buses_branches[bus_label] :
                     incident(g, bus_label)
            for e in _edges
                bus2 = opposite(e, bus_label)
                _r_add_bus!(h, buses_to_rm, edges_to_rm, g, _buses_branches, bus2)
                h[e...] = g[e...]
                push!(edges_to_rm, e)
            end
        end
    end

    A = incidence_matrix(g; oriented=true)

    sub = _initgraph()
    buses_to_rm = VLabel[]
    edges_to_rm = ELabel[]
    _r_add_bus!(sub, buses_to_rm, edges_to_rm, g, buses_branches, first(keys(buses_branches)))
    remaining = copy(g)
    foreach(e -> delete!(remaining, e[1], e[2]), edges_to_rm)
    foreach(b -> delete!(remaining, b), buses_to_rm)
    sub, remaining
end

function extract_export_phases(g::MetaGraph, phases::Vector{Float64}, buses::Vector{String})
    [phases[code_for(g, bus)] for bus in buses]
end

function extract_export_phases(g::MetaGraph, phases::JuMP.Containers.DenseAxisArray{Float64}, buses::Vector{String})
    [phases[bus] for bus in buses]
end

function equivalent_parameters(
    orig::MetaGraph,
    label::String,
    buses::AbstractArray{VLabel},
    include_branches_limites::Bool;
    outages::AbstractArray{ELabel}=ELabel[],)

    nb_fixed_buses = length(buses)
    sloped_g = copy(orig)
    foreach(bus -> (sloped_g[bus] *= sloped_g[bus] ≥ 0 ? 1 : 2), labels(sloped_g))
    _, _, stressed_gen = dc_flow_optim(sloped_g; outages=outages, slack_buses=buses, slack_phases=zeros(Float64, nb_fixed_buses))
    flows, _, injections = dc_flow_optim(orig; outages=outages, slack_buses=buses, slack_phases=zeros(Float64, nb_fixed_buses))
    internal_branches = collect(eachindex(flows))
    f_max = [orig[e...].p_max for e in eachindex(flows)]
    gen_slope = stressed_gen.data .- injections.data
    B = zeros(Float64, nb_fixed_buses, nb_fixed_buses)
    F = zeros(Float64, ne(orig), nb_fixed_buses)
    for i in 1:nb_fixed_buses
        dϕ = zeros(Float64, nb_fixed_buses)
        dϕ[i] = 1
        flows2, _, injections2 = dc_flow_optim(orig; outages=outages, slack_buses=buses, slack_phases=dϕ)
        B[:, i] .= injections2.data .- injections.data
        F[:, i] .= [flows2[k] - flows[k] for k in eachindex(flows)]
    end
    orig_outages = copy(orig)
    foreach(e -> delete!(orig_outages, e...), outages)
    cc_orig = connected_components(orig_outages)
    ccs = OTS_cc[]
    for cc in cc_orig
        cc_buses = VLabel[bus for bus in buses if code_for(orig_outages, bus) in cc]
        l = sum([max(orig[bus], 0) for bus in cc_buses])
        g = sum([-min(orig[bus], 0) for bus in cc_buses])
        push!(ccs, OTS_cc(cc_buses, l, g))
    end
    if include_branches_limites
        return Equivalent(label, buses, injections.data, gen_slope, B, ccs, internal_branches, [flows[k] for k in eachindex(flows)], f_max, F)
    else
        return Equivalent(label, buses, injections.data, gen_slope, B, ccs)
    end
end
