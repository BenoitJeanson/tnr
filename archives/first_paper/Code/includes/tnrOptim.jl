using Logging
using JuMP

include("tnr.jl")

function init_model(g::MetaGraph,
    contingencies::AbstractArray{Int}=Int[],
    is_single_ρ::Bool=true,
    ρ_min_bound::Float64=0.0,
    n_1_connectedness::Bool=true,
    bus_confs::Vector{BusConf}=BusConf[],
    allow_branch_openings::Bool=true,
    OTS_only::Bool=false,
    tnr_pf::TNR_PF_TYPE=tnr_pf_phase)

    # model = direct_model(Gurobi.Optimizer())
    model = Model(Gurobi.Optimizer)
    # set_silent(model)
    set_optimizer_attribute(model, "DualReductions", 0)
    # set_optimizer_attribute(model, "LogFile", joinpath("_research", "tmp", "my_log_file.txt"))
    model.ext[:r] = TNR(g, contingencies, is_single_ρ, ρ_min_bound, n_1_connectedness, bus_confs, allow_branch_openings, OTS_only, tnr_pf)

    model, model.ext[:r]
end

function create_variables!(model::Model, r::TNR)

    @variable(model, c_w[cases(r), edge_ids(r)], Bin)
    @variable(model, flows[cases(r), edge_ids(r)])
    @variable(model, overload)

    @variable(model, gen[cases(r), 1:nv(g)] ≥ 0) #TODO : change 1:nv(g)...
    @variable(model, load[cases(r), 1:nv(g)] ≥ 0)
    @variable(model, lostload[n_1cases(r)])

    @variable(model, c_flows[(r.n_1_connectedness ? 1 : nb_cases(r)):nb_cases(r), edge_ids(r)])

    @variable(model, cn1_a[n_1cases(r)], Bin)
    @variable(model, cn1_b[n_1cases(r)], Bin)#lower_bound = 0, upper_bound = 1)
    @variable(model, cn1_flows[n_1cases(r), edge_ids(r)])
    @variable(model, cn1_p_orig[n_1cases(r)])
    @variable(model, cn1_π[n_1cases(r), buses(r)], Bin)# lower_bound = 0, upper_bound = 1)
    @variable(model, cn1_ψ[n_1cases(r), bus in buses(r), incident(r, bus)], Bin)# lower_bound = 0, upper_bound = 1)

    if r.is_single_ρ
        @variable(model, ρ ≥ 0)
    else
        @variable(model, ρ[edge_ids(r)] ≥ 0)
    end

    if r.allow_branch_openings || r.OTS_only
        @variable(model, v_branch[1:ne(g)], Bin)
        if r.tnr_pf == tnr_pf_pst
            @variable(model, γ_branch[cases(r), edge_ids(r)])
        end
    end

    if r.tnr_pf == tnr_pf_pst
        big_M = 5 # TODO: Change
        @variable(model, -big_M .≤ β[n_1cases(r)] .≤ big_M)
    else
        @variable(model, ϕ[cases(r), buses(r)])
    end


    if !r.OTS_only
        @variable(model, hatgen[cases(r), bus_conf_ids(r)] ≥ 0)
        @variable(model, hatload[cases(r), bus_conf_ids(r)] ≥ 0)
        @variable(model, v_bus[bus_conf_ids(r)], Bin)
        @variable(model, γ_bus[cases(r), k=bus_conf_ids(r)])

        @variable(model, cn1_χ[n_1cases(r), bus in buses(r), incident(r, bus)], Bin)
        @variable(model, cn1_U[bus in keys(r.bus_to_conf_ids), branch in incident(r, bus); (bus, branch) in keys(r.bus_branch_to_conf_ids)], Bin)
        @variable(model, cn1_χ_circ[n_1cases(r), bus in buses(r), incident(r, bus)], Bin)
        @variable(model, cn1_𝚿[n_1cases(r), bus in buses(r), incident(r, bus)], Bin)
        @variable(model, cn1_hat_𝚿[n_1cases(r), bus_conf_ids(r)], Bin)
    end
end

function OTS_flows_phases!(model, r::TNR)
    gen = model[:gen]
    load = model[:load]
    flows = model[:flows]
    ϕ = model[:ϕ]
    c_w = model[:c_w]

    big_M = 10 #TODO: to change    bus_orig = 12#TODO: to change

    @constraint(model, [c in cases(r)], ϕ[c, bus_orig] == 0)
    @constraint(model, [c in cases(r), bus in buses(r); bus ≠ bus_orig], #TODO: probably bus ≠ bus_orig is not necessary but can help in case of numerical instability: the slack takes it.
        load[c, bus] - gen[c, bus] == sum(r.A'[bus, e] * flows[c, e] for e in edge_ids(r)))

    @constraint(model, [c in cases(r)], flows[c, :] .≤ (1 .- c_w[c, :]) .* big_M)
    @constraint(model, [c in cases(r)], -flows[c, :] .≤ (1 .- c_w[c, :]) .* big_M)

    B = [br.b for br in r.branches]

    @constraint(model, [c in cases(r), edge in edge_ids(r)],
        flows[c, edge] - sum(B[edge] * r.A'[bus, edge] * ϕ[c, bus] for bus in buses(r)) ≤ big_M * c_w[c, edge])
    @constraint(model, [c in cases(r), edge in edge_ids(r)],
        flows[c, edge] - sum(B[edge] * r.A'[bus, edge] * ϕ[c, bus] for bus in buses(r)) ≥ -big_M * c_w[c, edge])
end

function OTS_flows!(model, r::TNR)
    gen = model[:gen]
    load = model[:load]
    β = model[:β]
    flows = model[:flows]
    c_w = model[:c_w]

    big_M = 5 # TODO: change

    γ_branch = model[:γ_branch]

    @constraint(model, [c in cases(r)], γ_branch[c, :] .≤ c_w[c, :] .* big_M)
    @constraint(model, [c in cases(r)], -γ_branch[c, :] .≤ c_w[c, :] .* big_M)

    @constraint(model, [c in cases(r)], flows[c, :] .≤ (1 .- c_w[c, :]) .* big_M)
    @constraint(model, [c in cases(r)], -flows[c, :] .≤ (1 .- c_w[c, :]) .* big_M)

    @constraint(model, [c in n_1cases(r)], flows[c, r.outages[c]] == 0)  # nb_cases is the index of the base case

    D = spdiagm([e_index_for(r.g, e).b for e in edges(r.g)])
    B = r.A' * D * r.A
    N = nv(r.g)
    B_inv = inv(Matrix(B + fill(1 / N, N, N))) - fill(1 / N, N, N)
    p_to_f = D * r.A * B_inv
    δ_to_f = D - p_to_f * r.A' * D

    for c in cases(r)
        outage_phase = @expression(model, c == nb_cases(r) ? 0 : β[c] .* [is_outage(r, c, e) for e in edge_ids(r)])
        @constraint(model, flows[c, :] .== p_to_f * (load[c, :] .- gen[c, :]) .+ δ_to_f * (outage_phase .+ γ_branch[c, :]))
    end
end

function TNR_flows!(model, r::TNR)
    r.allow_branch_openings && (v_branch = model[:v_branch])
    gen = model[:gen]
    load = model[:load]
    hatgen = model[:hatgen]
    hatload = model[:hatload]
    v_bus = model[:v_bus]
    γ_bus = model[:γ_bus]
    β = model[:β]
    flows = model[:flows]

    big_M = 5 # TODO: change

    if r.allow_branch_openings
        γ_branch = model[:γ_branch]
        v_branch = model[:v_branch]

        @constraint(model, [c in cases(r)], γ_branch[c, :] .≤ v_branch .* big_M)
        @constraint(model, [c in cases(r)], -γ_branch[c, :] .≤ v_branch .* big_M)

        @constraint(model, [c in cases(r)], flows[c, :] .≤ (1 .- v_branch) .* big_M)
        @constraint(model, [c in cases(r)], -flows[c, :] .≤ (1 .- v_branch) .* big_M)
    end

    @constraint(model, [c in cases(r), bc_id in bus_conf_ids(r)], γ_bus[c, bc_id] ≤ v_bus[bc_id] * big_M)
    @constraint(model, [c in cases(r), bc_id in bus_conf_ids(r)], -γ_bus[c, bc_id] ≤ v_bus[bc_id] * big_M)

    @constraint(model, [c in n_1cases(r)], flows[c, r.outages[c]] == 0)  # nb_cases is the index of the base case

    D = spdiagm([e_index_for(r.g, e).b for e in edges(r.g)])
    B = r.A' * D * r.A
    N = nv(r.g)
    B_inv = inv(Matrix(B + fill(1 / N, N, N))) - fill(1 / N, N, N)
    p_to_f = D * r.A * B_inv
    δ_to_f = D - p_to_f * r.A' * D

    for c in cases(r)
        outage_phase = @expression(model, c == nb_cases(r) ? 0 : β[c] .* [is_outage(r, c, e) for e in edge_ids(r)])
        line_opening_phase = @expression(model, r.allow_branch_openings ? γ_branch[c, :] : zeros(edge_ids(r)))
        bus_splitting_phase = @expression(model, isempty(r.bus_confs) ? 0
                                                 : [sum(br in sb.branch_ids ? (r.A[br, bc.bus] * γ_bus[c, k]) : 0
                                                      for (k, bc) in enumerate(r.bus_confs) for sb in r.bus_confs[k].subBuses) for br = edge_ids(r)])

        @constraint(model, flows[c, :] .== p_to_f * (load[c, :] .- gen[c, :])
                                           +
                                           δ_to_f * (outage_phase .+ line_opening_phase .+ bus_splitting_phase))
    end

    for (bc_id, bus, sb) in ((bc_id, bc.bus, sb) for (bc_id, bc) in enumerate(r.bus_confs) for sb in bc.subBuses)
        @constraint(model, [c in cases(r)],
            sum(flows[c, br] * r.A[br, bus] for br in sb.branch_ids) - (hatload[c, bc_id] - hatgen[c, bc_id]) ≤ (1 - v_bus[bc_id]) * big_M)
        @constraint(model, [c in cases(r)],
            -(sum(flows[c, br] * r.A[br, bus] for br in sb.branch_ids) - (hatload[c, bc_id] - hatgen[c, bc_id])) ≤ (1 - v_bus[bc_id]) * big_M)
    end

end

"""
overload is the value to be minimized. It relies on ρ which is flow/flow_limit - ρ_min_bound. ρ_min_bound=1 if the expectation is to get a situation without overload. It can be either the max of load over all the branches (is_single_ρ) or it can be the sum of the of all the over the min_bound load.
"""
function overload!(model, r::TNR)
    flows = model[:flows]
    overload = model[:overload]
    ρ = model[:ρ]

    p_max = [e_index_for(r.g, br).p_max for br in edges(r.g)]
    if r.is_single_ρ
        @constraint(model, overload == ρ)

        @constraint(model, [i in cases(r)], flows[i, :] .≤ (ρ + r.ρ_min_bound) .* p_max)
        @constraint(model, [i in cases(r)], flows[i, :] .≥ -(ρ + r.ρ_min_bound) .* p_max)
    else
        @constraint(model, overload == sum(ρ))

        @constraint(model, [i in cases(r)], flows[i, :] .≤ (ρ .+ r.ρ_min_bound) .* p_max)
        @constraint(model, [i in cases(r)], flows[i, :] .≥ -(ρ .+ r.ρ_min_bound) .* p_max)
    end
end

function simple_overload!(model, r::TNR)
    flows = model[:flows]
    p_max = [e_index_for(r.g, br).p_max for br in edges(r.g)]
    @constraint(model, [i in cases(r)], flows[i, :] .≤ p_max)
    @constraint(model, [i in cases(r)], -flows[i, :] .≤ p_max)
end

function c_w!(model, r::TNR)
    c_w = model[:c_w]
    if r.allow_branch_openings || r.OTS_only
        v_branch = model[:v_branch]
        @constraint(model, [c in cases(r), e in edge_ids(r)], c_w[c, e] ≥ is_outage(r, c, e))
        @constraint(model, [c in cases(r), e in edge_ids(r)], c_w[c, e] ≥ v_branch[e])
        @constraint(model, [c in cases(r), e in edge_ids(r)], c_w[c, e] ≤ is_outage(r, c, e) + v_branch[e])
    else
        @constraint(model, [c in cases(r), e in edge_ids(r)], c_w[c, e] == is_outage(r, c, e))
    end
end

function OTS_N_connectedness!(model, r::TNR)
    c_w = model[:c_w]
    c_flows = model[:c_flows]

    big_M2 = nb_buses(r)
    cases = (r.n_1_connectedness ? 1 : nb_cases(r)):nb_cases(r)        # nb_cases is the index of the base case

    bus = (e, s) -> findfirst(v -> r.A[e, v] == s, buses(r))

    # # c_flows are null when a line is open : c_w == 1
    @constraint(model, [c in cases, e in edge_ids(r)], c_flows[c, e] ≤ big_M2 * (1 - c_w[c, e]))
    @constraint(model, [c in cases, e in edge_ids(r)], -c_flows[c, e] ≤ big_M2 * (1 - c_w[c, e]))

    # # the first bus holds a generator with n_bus_tot - 1 - its number of possible confs,
    # # each other consumes 1 + (the number of possible confs on the bus)
    # # This ensures the c_flows are balanced at the substation level
    @constraint(model, [c in cases, bus in buses(r)],
        (bus == 1 ? -nb_buses(r) : 0) + 1 - sum(c_flows[c, br] * r.A[br, bus] for br in edge_ids(r)) == 0)
end

function TNR_N_connectedness!(model, r::TNR)
    v_bus = model[:v_bus]
    c_w = model[:c_w]
    c_flows = model[:c_flows]

    nb_bus_tot = nb_buses(r) + (isempty(r.bus_confs) ? 0 : sum(length(r.bus_to_conf_ids[bus]) for bus in keys(r.bus_to_conf_ids)))
    big_M2 = nb_bus_tot + 1
    cases = (r.n_1_connectedness ? 1 : nb_cases(r)):nb_cases(r)        # nb_cases is the index of the base case

    # # c_flows are null when a line is open : c_w == 1
    @constraint(model, [c in cases, e in edge_ids(r)], c_flows[c, e] ≤ big_M2 * (1 - c_w[c, e]))
    @constraint(model, [c in cases, e in edge_ids(r)], -c_flows[c, e] ≤ big_M2 * (1 - c_w[c, e]))

    # # the first bus holds a generator with n_bus_tot - 1 - its number of possible confs,
    # # each other consumes 1 + (the number of possible confs on the bus)
    # # This ensures the c_flows are balanced at the substation level
    nb_conf = bus -> bus_nb_confs(r, bus)
    @constraint(model, [c in cases, bus in buses(r)],
        (bus == 1 ? -nb_bus_tot : 0) + 1 + nb_conf(bus) - sum(c_flows[c, br] * r.A[br, bus] for br in edge_ids(r)) == 0)

    # But when confs are activated, the balance needs also to be ensured for the concerned subBuses
    for (bc_id, bc) in enumerate(r.bus_confs), c in cases
        c_flows_on_bus_bus = @expression(model,
            sum(c_flows[c, br] * r.A[br, bc.bus] for sb in bc.subBuses for br in sb.branch_ids))
        @constraint(model, 1 - c_flows_on_bus_bus ≤ big_M2 * (1 - v_bus[bc_id]))
        @constraint(model, -(1 - c_flows_on_bus_bus) ≤ big_M2 * (1 - v_bus[bc_id]))
    end
end

function OTS_N_1_connectednes!(model, r::TNR, bus_origin)
    bigM_nb_v = nv(g) + 2
    bigM_π = 2

    c_w = model[:c_w]
    cn1_a = model[:cn1_a]
    cn1_b = model[:cn1_b]
    cn1_flows = model[:cn1_flows]
    cn1_p_orig = model[:cn1_p_orig]
    cn1_π = model[:cn1_π]
    cn1_ψ = model[:cn1_ψ]

    # @constraint(model, [n_1cases(r)], cn1_a[:] .+ cn1_b[:] .== 1)

    @constraint(model, [c in n_1cases(r), e in edge_ids(r)], cn1_flows[c, e] ≤ bigM_nb_v * (1 - c_w[c, e]))
    @constraint(model, [c in n_1cases(r), e in edge_ids(r)], -cn1_flows[c, e] ≤ bigM_nb_v * (1 - c_w[c, e]))

    @constraint(model, [c in n_1cases(r), bus in [bus for bus in buses(r) if bus ≠ bus_origin]],
        cn1_π[c, bus] - sum(r.A[e, bus] * cn1_flows[c, e] for e in edge_ids(r)) == 0)
    @constraint(model, [c in n_1cases(r)],
        cn1_p_orig[c] - sum(r.A[e, bus_origin] * cn1_flows[c, e] for e in edge_ids(r)) == 0)

    # @constraint(model, [c in n_1cases(r)],   cn1_p_orig[c] + nb_buses(r) - 1  ≤ bigM_nb_v * cn1_a[c])
    # @constraint(model, [c in n_1cases(r)], -(cn1_p_orig[c] + nb_buses(r) - 1) ≤ bigM_nb_v * cn1_a[c])

    @constraint(model, [c in n_1cases(r)], cn1_π[c, bus_origin] == 1)

    # @constraint(model, [c in n_1cases(r)],
    #       sum(abs(r.A[r.outages[c], bus]) * cn1_π[c, bus] for bus in buses(r)) - 1  ≤ bigM_π * cn1_b[c])
    # @constraint(model, [c in n_1cases(r)],
    #     -(sum(abs(r.A[r.outages[c], bus]) * cn1_π[c, bus] for bus in buses(r)) - 1) ≤ bigM_π * cn1_b[c])

    @constraint(model, [c in n_1cases(r), bus in buses(r); bus ≠ bus_origin],
        cn1_π[c, bus] ≤ sum(cn1_ψ[c, bus, e] for e in incident(r, bus)))

    @constraint(model, [c in n_1cases(r), bus in buses(r), e in incident(r, bus)],
        cn1_π[c, bus] ≥ cn1_ψ[c, bus, e])

    @constraint(model, [c in n_1cases(r), bus in buses(r), e in incident(r, bus)],
        cn1_ψ[c, bus, e] ≥ cn1_π[c, opposite(r, e, bus)] - c_w[c, e])
    @constraint(model, [c in n_1cases(r), bus in buses(r), e in incident(r, bus)],
        cn1_ψ[c, bus, e] ≤ cn1_π[c, opposite(r, e, bus)])
    @constraint(model, [c in n_1cases(r), bus in buses(r), e in incident(r, bus)],
        cn1_ψ[c, bus, e] ≤ 1 - c_w[c, e])
end

function TNR_N_1_connectednes!(model, r::TNR, bus_origin)
    buses_wo_O = (k for k in buses(r) if k ≠ bus_origin)
    bigM_nb_v = nb_buses(r) + 2 + nb_bus_confs(r)

    v_bus = model[:v_bus]
    c_w = model[:c_w]
    cn1_a = model[:cn1_a]
    cn1_b = model[:cn1_b]
    cn1_flows = model[:cn1_flows]
    cn1_p_orig = model[:cn1_p_orig]
    cn1_π = model[:cn1_π]
    cn1_ψ = model[:cn1_ψ]

    cn1_U = model[:cn1_U]

    cn1_χ = model[:cn1_χ]
    cn1_χ_circ = model[:cn1_χ_circ]
    cn1_𝚿 = model[:cn1_𝚿]
    cn1_hat_𝚿 = model[:cn1_hat_𝚿]

    terminal = e -> [bus for bus in buses(r) if r.A[e, bus] ≠ 0]

    @variable(model, cn1_hatπ[n_1cases(r), bc in bus_conf_ids(r)], Bin)# lower_bound = 0, upper_bound = 1)

    for c in n_1cases(r), bus in buses(r), edge in incident(r, bus)
        @constraint(model, cn1_χ[c, bus, edge] == cn1_χ_circ[c, bus, edge]
                                                  +
                                                  sum(cn1_hatπ[c, conf]
                                                      for conf in get(r.bus_to_conf_ids, opposite(r, edge, bus), Int[])
                                                      if edge in r.bus_confs[conf].subBuses[1].branch_ids))
        @constraint(model, cn1_χ_circ[c, bus, edge] ≤ cn1_π[c, opposite(r, edge, bus)])
        @constraint(model, cn1_χ_circ[c, bus, edge] ≤ 1
                                                      -
                                                      sum(v_bus[conf]
                                                          for conf in get(r.bus_to_conf_ids, opposite(r, edge, bus), Int[])
                                                          if edge in r.bus_confs[conf].subBuses[1].branch_ids))
        @constraint(model, cn1_χ_circ[c, bus, edge] ≥ cn1_π[c, opposite(r, edge, bus)]
                                                      -
                                                      sum(v_bus[conf]
                                                          for conf in get(r.bus_to_conf_ids, opposite(r, edge, bus), Int[])
                                                          if edge in r.bus_confs[conf].subBuses[1].branch_ids))

        @constraint(model, cn1_ψ[c, bus, edge] ≤ cn1_χ[c, bus, edge])
        @constraint(model, cn1_ψ[c, bus, edge] ≤ 1 - c_w[c, edge])
        @constraint(model, cn1_ψ[c, bus, edge] ≥ cn1_χ[c, bus, edge] - c_w[c, edge])
    end

    @constraint(model, [c in n_1cases(r),
            bus in buses(r),
            conf in get(r.bus_to_conf_ids, bus, Int[]),
            edge in r.bus_confs[conf].subBuses[1].branch_ids],
        cn1_hat_𝚿[c, conf] ≥ cn1_ψ[c, bus, edge])
    @constraint(model, [c in n_1cases(r),
            bus in buses(r),
            conf in get(r.bus_to_conf_ids, bus, Int[])],
        cn1_hat_𝚿[c, conf] ≤
        sum(cn1_ψ[c, bus, edge] for edge in r.bus_confs[conf].subBuses[1].branch_ids))

    @constraint(model, [c in n_1cases(r), conf in bus_conf_ids(r)],
        cn1_hatπ[c, conf] ≤ v_bus[conf])
    @constraint(model, [c in n_1cases(r), conf in bus_conf_ids(r)],
        cn1_hatπ[c, conf] ≤ cn1_hat_𝚿[c, conf])
    @constraint(model, [c in n_1cases(r), conf in bus_conf_ids(r)],
        cn1_hatπ[c, conf] ≥ v_bus[conf] + cn1_hat_𝚿[c, conf] - 1)

    for ((bus, branch), confs) in r.bus_branch_to_conf_ids
        @constraint(model, [conf in confs], cn1_U[bus, branch] ≥ v_bus[conf])
        @constraint(model, cn1_U[bus, branch] ≤ sum(v_bus[conf] for conf in confs))
    end

    for bus in buses(r), edge in incident(r, bus)
        if (bus, edge) in keys(r.bus_branch_to_conf_ids)
            @constraint(model, [c in n_1cases(r)], cn1_𝚿[c, bus, edge] ≤ 1 - cn1_U[bus, edge])
            @constraint(model, [c in n_1cases(r)], cn1_𝚿[c, bus, edge] ≤ cn1_ψ[c, bus, edge])
            @constraint(model, [c in n_1cases(r)], cn1_𝚿[c, bus, edge] ≥ cn1_ψ[c, bus, edge] - cn1_U[bus, edge])
        else
            @constraint(model, [c in n_1cases(r)], cn1_𝚿[c, bus, edge] == cn1_ψ[c, bus, edge])
        end
    end

    @constraint(model, [c in n_1cases(r), bus in buses(r), edge in incident(r, bus)],
        cn1_π[c, bus] ≥ cn1_𝚿[c, bus, edge])
    @constraint(model, [c in n_1cases(r), bus in buses(r)],
        cn1_π[c, bus] ≤ sum(cn1_𝚿[c, bus, edge] for edge in incident(r, bus)))

    @constraint(model, cn1_a[:] .+ cn1_b[:] .== 1)

    @constraint(model, [c in n_1cases(r)], sum(cn1_χ[c, term, r.outages[c]] for term in terminal(r.outages[c])) - 1 ≤ 1 - cn1_b[c])
    @constraint(model, [c in n_1cases(r)], -(sum(cn1_χ[c, term, r.outages[c]] for term in terminal(r.outages[c])) - 1) ≤ 1 - cn1_b[c])

    n_1_sum_us = @expression(model, nv(g) - 1 + sum(v_bus[conf] for conf in bus_conf_ids(r)))
    @constraint(model, [n_1cases(r)], cn1_p_orig[:] .- n_1_sum_us .≤ bigM_nb_v .* (1 .- cn1_a[:]))
    @constraint(model, [n_1cases(r)], cn1_p_orig[:] .- n_1_sum_us .≥ -bigM_nb_v .* (1 .- cn1_a[:]))

    @constraint(model, [c in n_1cases(r), e in edge_ids(r)], cn1_flows[c, e] ≤ bigM_nb_v * (1 - c_w[c, e]))
    @constraint(model, [c in n_1cases(r), e in edge_ids(r)], -cn1_flows[c, e] ≤ bigM_nb_v * (1 - c_w[c, e]))

    @constraint(model, [c in n_1cases(r), bus_conf_ids(r)], cn1_hatπ[c, :] .≤ v_bus[:])

    for (bc_id, bus, sb) in ((bc_id, bc.bus, sb) for (bc_id, bc) in enumerate(r.bus_confs) for sb in bc.subBuses)
        @constraint(model, [c in n_1cases(r)],
            sum(cn1_flows[c, edge] * r.A[edge, bus] for edge in sb.branch_ids) - cn1_hatπ[c, bc_id] ≤ (1 - v_bus[bc_id]) * bigM_nb_v)
        @constraint(model, [c in n_1cases(r)],
            sum(cn1_flows[c, edge] * r.A[edge, bus] for edge in sb.branch_ids) - cn1_hatπ[c, bc_id] ≥ -(1 - v_bus[bc_id]) * bigM_nb_v)
    end

    @constraint(model, [c in n_1cases(r), bus in buses_wo_O],
        sum(cn1_hatπ[c, bc_id] for bc_id in get(r.bus_to_conf_ids, bus, Int[])) + cn1_π[c, bus] - sum(r.A[edge, bus] * cn1_flows[c, edge] for edge in incident(r, bus)) == 0)

    @constraint(model, [c in n_1cases(r)], cn1_π[c, bus_origin] == 1)

    @constraint(model, [c in n_1cases(r)],
        cn1_p_orig[c] + sum(r.A[edge, bus_origin] * cn1_flows[c, edge] for edge in incident(r, bus_origin)) == 0)
end

function N_balance!(model, r::TNR)
    load = model[:load]
    gen = model[:gen]

    for bus in buses(r)
        if r.p[bus] < 0
            @constraint(model, load[nb_cases(r), bus] == 0)
            @constraint(model, gen[nb_cases(r), bus] == -r.p[bus])
        else
            @constraint(model, gen[nb_cases(r), bus] == 0)
            @constraint(model, load[nb_cases(r), bus] == r.p[bus])
        end
    end
end

function OTS_balance!(model, r::TNR)
    bigM = 10 # TODO: to define

    gen = model[:gen]
    load = model[:load]

    @variable(model, σ[n_1cases(r)] ≥ 0)

    if haskey(model, :cn1_π)
        cn1_π = model[:cn1_π]
        for c in n_1cases(r), bus in buses(r)
            if r.p[bus] < 0
                @constraint(model, load[c, bus] == 0)
                @constraint(model, gen[c, bus] + σ[c] * r.p[bus] ≤ bigM * (1 - cn1_π[c, bus]))
                @constraint(model, -(gen[c, bus] + σ[c] * r.p[bus]) ≤ bigM * (1 - cn1_π[c, bus]))
                @constraint(model, gen[c, bus] ≤ bigM * cn1_π[c, bus])
                @constraint(model, -gen[c, bus] ≤ bigM * cn1_π[c, bus])
            else
                @constraint(model, gen[c, bus] == 0)
                @constraint(model, load[c, bus] - r.p[bus] ≤ bigM * (1 - cn1_π[c, bus]))
                @constraint(model, -(load[c, bus] - r.p[bus]) ≤ bigM * (1 - cn1_π[c, bus]))
                @constraint(model, load[c, bus] ≤ bigM * cn1_π[c, bus])
                @constraint(model, -load[c, bus] ≤ bigM * cn1_π[c, bus])
            end
        end
    end

    @constraint(model, [c in n_1cases(r)], sum(load[c, :] .- gen[c, :]) == 0)
end

function TNR_balance!(model, r::TNR)
    bigM = 2 * maximum(abs, r.p)

    v_bus = model[:v_bus]
    cn1_π = model[:cn1_π]
    cn1_hatπ = model[:cn1_hatπ]

    gen = model[:gen]
    load = model[:load]
    hatgen = model[:hatgen]
    hatload = model[:hatload]

    @variable(model, σ[n_1cases(r)] ≥ 0)

    @variable(model, cn1_hatd_u[bus_conf_ids(r)] ≥ 0)
    @variable(model, cn1_hatg_u[n_1cases(r), bus_conf_ids(r)] ≥ 0)

    hatp = [bc.subBuses[1].p for bc in r.bus_confs]

    # N connectedness is supposed to be ensured. hat_p_s = u_s * p_conf_s
    for bc_id in bus_conf_ids(r)
        _p = hatp[bc_id]
        if _p > 0
            @constraint(model, hatload[nb_cases(r), bc_id] ≤ bigM * v_bus[bc_id])
            @constraint(model, hatload[nb_cases(r), bc_id] - _p ≤ bigM * (1 - v_bus[bc_id]))
            @constraint(model, -(hatload[nb_cases(r), bc_id] - _p) ≤ bigM * (1 - v_bus[bc_id]))
            @constraint(model, hatgen[nb_cases(r), bc_id] == 0)
        elseif _p < 0
            @constraint(model, -hatgen[nb_cases(r), bc_id] ≤ bigM * v_bus[bc_id])
            @constraint(model, hatgen[nb_cases(r), bc_id] + _p ≤ bigM * (1 - v_bus[bc_id]))
            @constraint(model, -(hatgen[nb_cases(r), bc_id] + _p) ≤ bigM * (1 - v_bus[bc_id]))
            @constraint(model, hatload[nb_cases(r), bc_id] == 0)
        else
            @constraint(model, hatgen[nb_cases(r), bc_id] == 0)
            @constraint(model, hatload[nb_cases(r), bc_id] == 0)
        end
    end

    for c in n_1cases(r)
        @constraint(model, sum(gen[c, bus] - load[c, bus] for bus in buses(r))
                           ==
                           0)
        for bc_id in bus_conf_ids(r)
            _p = hatp[bc_id]
            if _p > 0
                @constraint(model, hatload[c, bc_id] ≤ bigM * cn1_hatπ[c, bc_id])
                @constraint(model, hatload[c, bc_id] - _p ≤ bigM * (1 - cn1_hatπ[c, bc_id]))
                @constraint(model, -(hatload[c, bc_id] - _p) ≤ bigM * (1 - cn1_hatπ[c, bc_id]))
                @constraint(model, hatgen[c, bc_id] == 0)
            elseif _p < 0
                # println("c:$c, bc_id:$bc_id hatgen:$hatgen")
                @constraint(model, hatgen[c, bc_id] ≤ bigM * cn1_hatπ[c, bc_id])
                @constraint(model, hatgen[c, bc_id] + _p * σ[c] ≤ bigM * (1 - cn1_hatπ[c, bc_id]))
                @constraint(model, -(hatgen[c, bc_id] + _p * σ[c]) ≤ bigM * (1 - cn1_hatπ[c, bc_id]))
                @constraint(model, hatload[c, bc_id] == 0)
            else
                @constraint(model, hatgen[c, bc_id] == 0)
                @constraint(model, hatload[c, bc_id] == 0)
            end
        end
    end

    for bc_id in bus_conf_ids(r)
        _p = hatp[bc_id]
        if _p > 0
            @constraint(model, [c in n_1cases(r)], cn1_hatg_u[c, bc_id] == 0)
            @constraint(model, cn1_hatd_u[bc_id] ≤ bigM * v_bus[bc_id])
            @constraint(model, cn1_hatd_u[bc_id] - _p ≤ bigM * (1 - v_bus[bc_id]))
            @constraint(model, -(cn1_hatd_u[bc_id] - _p) ≤ bigM * (1 - v_bus[bc_id]))
        elseif _p < 0
            @constraint(model, cn1_hatd_u[bc_id] == 0)
            for c in n_1cases(r)
                @constraint(model, cn1_hatg_u[c, bc_id] ≤ bigM * v_bus[bc_id])
                @constraint(model, cn1_hatg_u[c, bc_id] + σ[c] * _p ≤ bigM * (1 - v_bus[bc_id]))
                @constraint(model, -(cn1_hatg_u[c, bc_id] + σ[c] * _p) ≤ bigM * (1 - v_bus[bc_id]))
            end
        else
            @constraint(model, [c in n_1cases(r)], cn1_hatg_u[c, bc_id] == 0)
            @constraint(model, cn1_hatd_u[bc_id] == 0)
        end
    end

    for c in n_1cases(r), bus in buses(r)
        bc_ids = get(r.bus_to_conf_ids, bus, Int[])
        @constraint(model, load[c, bus] - sum(hatload[c, bc_id] for bc_id in bc_ids) ≤ bigM * cn1_π[c, bus])
        @constraint(model, -(load[c, bus] - sum(hatload[c, bc_id] for bc_id in bc_ids)) ≤ bigM * cn1_π[c, bus])
        @constraint(model, gen[c, bus] - sum(hatgen[c, bc_id] for bc_id in bc_ids) ≤ bigM * cn1_π[c, bus])
        @constraint(model, -(gen[c, bus] - sum(hatgen[c, bc_id] for bc_id in bc_ids)) ≤ bigM * cn1_π[c, bus])
        @constraint(model, load[c, bus] - maximum((0, r.p[bus]))
                           +
                           sum(cn1_hatd_u[bc_id] - hatload[c, bc_id] for bc_id in bc_ids) ≤ bigM * (1 - cn1_π[c, bus]))
        @constraint(model, -(load[c, bus] - maximum((0, r.p[bus]))
                             +
                             sum(cn1_hatd_u[bc_id] - hatload[c, bc_id] for bc_id in bc_ids)) ≤ bigM * (1 - cn1_π[c, bus]))
        @constraint(model, gen[c, bus] + σ[c] * minimum((0, r.p[bus]))
                           + sum(cn1_hatg_u[c, bc_id] - hatgen[c, bc_id] for bc_id in bc_ids) ≤ bigM * (1 - cn1_π[c, bus]))
        @constraint(model, -(gen[c, bus] + σ[c] * minimum((0, r.p[bus]))
                             + sum(cn1_hatg_u[c, bc_id] - hatgen[c, bc_id] for bc_id in bc_ids)) ≤ bigM * (1 - cn1_π[c, bus]))
    end
end

function loadloss(model, r::TNR)
    lostload = model[:lostload]
    load = model[:load]

    @constraint(model, [c in n_1cases(r)], lostload[c] == sum(r.p[bus] - load[c, bus] for bus in buses(r) if r.p[bus] ≥ 0))
end

function secured_dc_OTS(g::MetaGraph;
    contingencies::AbstractArray{Int}=Int[],
    is_single_ρ::Bool=true,
    ρ_min_bound::Float64=0.0,
    n_1_connectedness::Bool=false,
    bus_confs::Vector{BusConf}=BusConf[],
    allow_branch_openings::Bool=true,
    OTS_only::Bool=false,
    tnr_pf::TNR_PF_TYPE=tnr_pf_phase)

    model, r = init_model(g, contingencies, is_single_ρ, ρ_min_bound, n_1_connectedness, bus_confs, allow_branch_openings, OTS_only, tnr_pf)
    create_variables!(model, r)

    c_w!(model, r)
    N_balance!(model, r)

    if r.OTS_only
        if r.tnr_pf == tnr_pf_pst
            OTS_flows!(model, r)
        elseif r.tnr_pf == tnr_pf_phase
            OTS_flows_phases!(model, r)
        end
        OTS_N_connectedness!(model, r)

        OTS_N_1_connectednes!(model, r, 12)# TODO: improve bus_orig choice !

        OTS_balance!(model, r)
    else
        TNR_flows!(model, r)
        TNR_N_connectedness!(model, r)
        TNR_N_1_connectednes!(model, r, 1)
        TNR_balance!(model, r)

        @constraint(model, [bus in keys(r.bus_to_conf_ids)], sum(model[:v_bus][bc] for bc in r.bus_to_conf_ids[bus]) ≤ 1)
    end

    loadloss(model, r)

    simple_overload!(model, r)
    # overload!(model, r)
    # overload = model[:overload]
    # @constraint(model, overload ≤ 0)

    lostload = model[:lostload]
    @objective(model, Min,
        sum(lostload[c] for c in n_1cases(r)))
    # overload)
    # + 0.01  * (allow_branch_openings ? sum(model[:v_branch]) : 0)) 
    # + 0.01 * (!isempty(bus_confs)   ? sum(model[:v_bus])    : 0) )
    model, r
end