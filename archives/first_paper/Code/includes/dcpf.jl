using Logging
using SparseArrays
using MetaGraphsNext
using JuMP

@enum DCPF_Type pf_linalg pf_phase

function _incidence_matrix(g::MetaGraph; outages = [], trip = nothing)
    A = incidence_matrix(g; oriented = true)
    for outage in outages
        A[:, outage,] .= 0
    end
    if trip ≠ nothing
        A[:, trip] .= 0
    end
    return A
end


function _graph_to_mat(g::MetaGraph; outages = [], trip = nothing)
    A = _incidence_matrix(g, outages = outages, trip = trip)
    B = [e_index_for(g,e).b for e in edges(g)]
    P = [g[label_for(g, v)] for v in vertices(g)]
    A, B, P
end


function setflows!(g::MetaGraph, flows::Vector{Float64})
    for (i,e) in enumerate(edges(g))
        e_index_for(g, e).p = flows[i]
    end
end


function dc_flow(g::MetaGraph; outages = [], trip = nothing, pf_type::DCPF_Type = pf_linalg, kwargs...)
    h = deepcopy(g)
    dc_flow!(h, outages = outages, trip = trip, pf_type = pf_type)
end


function dc_flow!(g::MetaGraph; outages = [], trip = nothing, pf_type::DCPF_Type = pf_linalg, kwargs...)
    A, B, P = _graph_to_mat(g, outages = outages, trip = trip)
    if pf_type == pf_linalg
        flows = dc_flow_linalg(A, B, P)
    elseif pf_type == pf_phase
        flows = dc_flow_phase(A, B, P, kwargs...)
    end
    setflows!(g, flows)
    g
end


function dc_flow_linalg(A::AbstractMatrix, B::Vector{Float64}, P::Vector{Float64})
    @info "DC PF with linalg"
    _D = spdiagm(B)
    _B = A * _D * A'
    n = length(P)
    _B_inv = inv(Matrix(_B + fill(1/n, n, n))) - fill(1/n, n, n)
    p_to_f = _D * A' * _B_inv
    flows = p_to_f * P
end


function dc_flow_phase(A::AbstractMatrix, B::Vector{Float64}, P::Vector{Float64}; kwargs...)
    @info "DC PF with optim phase"
    nb_buses, nb_edges = size(A)
    slack = get(kwargs, :slack, findmin(P)[2])
    model = Model(Gurobi.Optimizer);
    set_silent(model)
    @variable(model, ϕ[1:nb_buses])
    @variable(model, flows[1:nb_edges])
    @constraint(model, ϕ[slack] == 0)
    @constraint(model, [bus in 1:nb_buses; bus ≠ slack], P[bus] == sum(A[bus, e] * flows[e] for e in 1:nb_edges))
    @constraint(model, [edge in 1:nb_edges], flows[edge] == sum(B[edge] * A[bus, edge] * ϕ[bus] for bus in 1:nb_buses))
    optimize!(model);
    if is_solved_and_feasible(model)
        @info "dc_flow_phase solved"
        return value.(model[:flows])
    else
        @warn "dc_flow_phase is not solved"
        return zeros(nb_edges)
    end
end