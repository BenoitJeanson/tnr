using Logging
using SparseArrays
using MetaGraphsNext
using JuMP
using Gurobi

include("../src/GraphUtils.jl")
include("../wip/Equivalent.jl")

@enum DCPF_Type pf_linalg pf_optim

function _incidence_matrix(g::MetaGraph; outages=[], trip=nothing)
    A = incidence_matrix(g; oriented=true)
    for outage in outages
        A[:, outage,] .= 0
    end
    if trip ≠ nothing
        A[:, trip] .= 0
    end
    return A
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

function setflows!(g::MetaGraph, flows::JuMP.Containers.SparseAxisArray{Float64,2,ELabel})
    for edge in eachindex(flows)
        g[edge...].p = flows[edge]
    end
end


function dc_flow(g::MetaGraph; outages=[], trip=nothing, pf_type::DCPF_Type=pf_optim, kwargs...)
    h = deepcopy(g)
    dc_flow!(h; outages=outages, trip=trip, pf_type=pf_type, kwargs...)
end


function dc_flow!(g::MetaGraph; outages=[], trip=nothing, pf_type::DCPF_Type=pf_optim, kwargs...)
    A, B, P = _graph_to_mat(g, outages=outages, trip=trip)
    if pf_type == pf_linalg
        flows = dc_flow_linalg(A, B, P)
    elseif pf_type == pf_optim
        flows, _, _ = dc_flow_optim(A, B, P; kwargs...)
    end
    setflows!(g, flows)
    g
end


function dc_flow_linalg(A::AbstractMatrix, B::Vector{Float64}, P::Vector{Float64})
    _D = spdiagm(B)
    _B = A * _D * A'
    n = length(P)
    _B_inv = inv(Matrix(_B + fill(1 / n, n, n))) - fill(1 / n, n, n)
    p_to_f = _D * A' * _B_inv
    flows = p_to_f * P
end


function dc_flow_optim(A::AbstractMatrix, B::Vector{Float64}, P::Vector{Float64};
    slack::Union{Nothing,Int}=nothing,
    slack_buses::Union{Nothing,Vector{Int}}=nothing,
    slack_phases::Union{Nothing,Vector{Float64}}=nothing,
    equivalent::Union{Nothing,Equivalent}=nothing)

    nb_buses, nb_edges = size(A)
    model = Model(Gurobi.Optimizer)
    set_silent(model)
    @variable(model, ϕ[1:nb_buses])
    @variable(model, flows[1:nb_edges])

    (_slack_buses, _slack_phases) = !isnothing(slack_buses) && !isnothing(slack_phases) ?
                                    (slack_buses, slack_phases) :
                                    ([isnothing(slack) ? findmin(P)[2] : slack], [0.0])

    @constraint(model, [f_bus_id in 1:length(_slack_buses)], ϕ[_slack_buses[f_bus_id]] == _slack_phases[f_bus_id])
    @constraint(model, [edge in 1:nb_edges], flows[edge] == sum(B[edge] * A[bus, edge] * ϕ[bus] for bus in 1:nb_buses))

    if isnothing(equivalent)
        @constraint(model, [bus in 1:nb_buses; !(bus in _slack_buses)],
            P[bus] ==
            sum(A[bus, e] * flows[e] for e in 1:nb_edges))
    else
        @variable(model, flows_e[equivalent.buses])
        @constraint(model, [bus in 1:nb_buses; !(bus in _slack_buses)],
            P[bus] ==
            sum(A[bus, e] * flows[e] for e in 1:nb_edges) -
            (bus in equivalent.buses ? flows_e[bus] : 0))
        @constraint(model, [e_bus_id in 1:length(equivalent.buses)],
            flows_e[equivalent.buses[e_bus_id]] ==
            equivalent.p[e_bus_id] +
            sum(equivalent.B[e_bus_id, other_id] * ϕ[other_bus] for (other_id, other_bus) in enumerate(equivalent.buses)))
        @constraint(model, [cc in equivalent.cc],
            sum(flows_e[bus] for bus in cc.buses) ==
            sum(equivalent.p[i] for (i, bus) in enumerate(equivalent.buses) if bus in cc.buses))
    end

    optimize!(model)

    if is_solved_and_feasible(model)
        slack_injections = [-sum(A[bus, e] * value(model[:flows][e]) for e in 1:nb_edges) for bus in _slack_buses]
        return value.(model[:flows]), value.(model[:ϕ]), slack_injections
    else
        @warn "dc_flow_phase is not solved"
        return zeros(nb_edges), zeros(nb_buses), isnothing(slack_phases) ? Int[] : zeros(nb_slack_ϕ)
    end
end

function dc_flow_optim(
    g::MetaGraph;
    outages::AbstractArray{ELabel}=ELabel[],
    trip=nothing,
    slack::Union{Nothing,VLabel}=nothing,
    slack_buses::Union{Nothing,Vector{VLabel}}=nothing,
    slack_phases::Union{Nothing,Vector{Float64}}=nothing,
    equivalent::Union{Nothing,Equivalent}=nothing)


    (_slack_buses, _slack_phases) = !isnothing(slack_buses) && !isnothing(slack_phases) ?
                                    (slack_buses, slack_phases) :
                                    ([isnothing(slack) ? label_for(g, findmin(bus -> g[bus], labels(g))[2]) : slack], [0.0])

    _outages = isnothing(trip) ? outages : [outages; trip]
    model = Model(Gurobi.Optimizer)
    set_silent(model)
    @variable(model, ϕ[collect(labels(g))])
    @variable(model, flows[bus1 in labels(g), bus2 in outneighbor_labels(g, bus1)])
    @variable(model, slack_injections[_slack_buses])

    @constraint(model, [bus in _slack_buses],
        slack_injections[bus] == sum(sign * flows[e...] for (e, sign) in incident_signed(g, bus)))

    @constraint(model, [(i, f_bus_id) in enumerate(_slack_buses)], ϕ[f_bus_id] == _slack_phases[i])
    @constraint(model, [edge in edge_labels(g); !(edge in _outages)], flows[edge...] == g[edge...].b * (ϕ[edge[2]] - ϕ[edge[1]]))
    @constraint(model, [edge in _outages], flows[edge...] == 0)

    if isnothing(equivalent)
        @constraint(model, [bus in labels(g); !(bus in _slack_buses)],
            g[bus] == sum(sign * flows[e...] for (e, sign) in incident_signed(g, bus)))
    else
        @variable(model, flows_e[equivalent.buses])
        @constraint(model, [bus in labels(g); !(bus in _slack_buses)],
            g[bus] ==
            sum(sign * flows[e...] for (e, sign) in incident_signed(g, bus)) +
            (bus in equivalent.buses ? flows_e[bus] : 0))
        @constraint(model, [e_bus_id in 1:length(equivalent.buses)],
            flows_e[equivalent.buses[e_bus_id]] ==
            equivalent.p[e_bus_id] +
            sum(equivalent.B[e_bus_id, other_id] * ϕ[other_bus] for (other_id, other_bus) in enumerate(equivalent.buses)))
        @constraint(model, [cc in equivalent.cc],
            sum(flows_e[bus] for bus in cc.buses) ==
            sum(equivalent.p[i] for (i, bus) in enumerate(equivalent.buses) if bus in cc.buses))
    end
    optimize!(model)

    if is_solved_and_feasible(model)
        return value.(model[:flows]), value.(model[:ϕ]), value.(model[:slack_injections])
    else
        @warn "dc_flow_phase is not solved"
        return nothing
    end
end


function dc_flow_optim!(
    g::MetaGraph;
    outages=ELabel[],
    trip=nothing,
    slack::Union{Nothing,VLabel}=nothing,
    slack_buses::Union{Nothing,Vector{VLabel}}=nothing,
    slack_phases::Union{Nothing,Vector{Float64}}=nothing,
    equivalent::Union{Nothing,Equivalent}=nothing)
    flows, ϕ, slack_injections = dc_flow_optim(g; outages, trip, slack, slack_buses, slack_phases, equivalent)
    setflows!(g, flows)
    flows, ϕ, slack_injections
end