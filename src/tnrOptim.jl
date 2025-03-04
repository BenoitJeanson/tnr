using Logging
using JuMP

include("tnr.jl")

function init_model(g::MetaGraph,
    contingencies::AbstractArray{Int},
    n_1_connectedness::Bool,
    bus_confs::Vector{BusConf},
    allow_branch_openings::Bool,
    OTS_only::Bool,
    tnr_pf::TNR_PF_TYPE,
    opf::Bool,
    bus_orig::String)

    # model = direct_model(Gurobi.Optimizer())
    model = Model(Gurobi.Optimizer)
    # set_silent(model)
    set_optimizer_attribute(model, "DualReductions", 0)
    # set_optimizer_attribute(model, "LogFile", joinpath("_research", "tmp", "my_log_file.txt"))
    model.ext[:r] = TNR(g, contingencies, n_1_connectedness, bus_confs, allow_branch_openings, OTS_only, tnr_pf, opf, bus_orig)

    model, model.ext[:r]
end

function create_variables!(m::Model, r::TNR)

    @variable(m, c_w[cases(r), edge_ids(r)], Bin)
    @variable(m, flows[cases(r), edge_ids(r)])
    @variable(m, overload)

    @variable(m, gen[cases(r), 1:nv(g)] ≥ 0) #TODO : change 1:nv(g)...
    @variable(m, load[cases(r), 1:nv(g)] ≥ 0)
    @variable(m, lostload[n_1cases(r)])

    @variable(m, c_flows[(r.n_1_connectedness ? 1 : nb_cases(r)):nb_cases(r), edge_ids(r)])

    @variable(m, cn1_flows[n_1cases(r), edge_ids(r)])
    @variable(m, cn1_p_orig[n_1cases(r)])
    @variable(m, cn1_π[n_1cases(r), buses(r)], Bin)# lower_bound = 0, upper_bound = 1)
    @variable(m, cn1_ψ[n_1cases(r), bus in buses(r), incident(r, bus)], Bin)# lower_bound = 0, upper_bound = 1)

    if r.allow_branch_openings || r.OTS_only
        @variable(m, v_branch[edge_ids(r)], Bin)
        if r.tnr_pf == tnr_pf_pst
            @variable(m, γ_branch[cases(r), edge_ids(r)])
        end
    end

    if r.tnr_pf == tnr_pf_pst
        big_M = 5 # TODO: Change
        @info "HARDCODED create_variables big_M for tnr-pf β boundaries: $big_M"
        @variable(m, -big_M .≤ β[n_1cases(r)] .≤ big_M)
    else
        @variable(m, ϕ[cases(r), buses(r)])
    end

    if !r.OTS_only
        @variable(m, hatgen[cases(r), bus_conf_ids(r)] ≥ 0)
        @variable(m, hatload[cases(r), bus_conf_ids(r)] ≥ 0)
        @variable(m, v_bus[bus_conf_ids(r)], Bin)

        @variable(m, cn1_hatπ[n_1cases(r), bc in bus_conf_ids(r)], Bin)# lower_bound = 0, upper_bound = 1)
        @variable(m, cn1_χ[n_1cases(r), bus in buses(r), incident(r, bus)], Bin)
        @variable(m, cn1_U[bus in keys(r.bus_to_conf_ids), branch in incident(r, bus); (bus, branch) in keys(r.bus_branch_to_conf_ids)], Bin)
        @variable(m, cn1_χ_circ[n_1cases(r), bus in buses(r), incident(r, bus)], Bin)
        @variable(m, cn1_𝚿[n_1cases(r), bus in buses(r), incident(r, bus)], Bin)
        @variable(m, cn1_hat_𝚿[n_1cases(r), bus_conf_ids(r)], Bin)

        if r.tnr_pf == tnr_pf_pst
            @variable(m, γ_bus[cases(r), k=bus_conf_ids(r)])
        elseif r.tnr_pf == tnr_pf_phase
            @variable(m, ϕ_e[cases(r), edge_ids(r), 1:2])
            @variable(m, ϕ_bc[cases(r), bus_conf_ids(r)])
        end
    end

    if r.opf
        @variable(m, loadshed[buses(r)] ≥ 0)
        @variable(m, σ[cases(r)] ≥ 0)
    else
        @variable(m, σ[n_1cases(r)] ≥ 0)
    end

end

function OTS_flows_phases!(m::Model, r::TNR)
    big_M = 1000 #TODO: to change  
    @info "HARDCODED OTS_flows_phases bigM to $big_M"

    @constraint(m, [c in cases(r)], m[:ϕ][c, r.bus_orig_id] == 0)
    @constraint(m, [c in cases(r), bus in buses(r); bus ≠ r.bus_orig_id], #TODO: probably bus ≠ bus_orig is not necessary but can help in case of numerical instability: the slack takes it.
        m[:load][c, bus] - m[:gen][c, bus] == sum(r.A'[bus, e] * m[:flows][c, e] for e in edge_ids(r)))

    @constraint(m, [c in cases(r)], m[:flows][c, :] .≤ (1 .- m[:c_w][c, :]) .* big_M)
    @constraint(m, [c in cases(r)], -m[:flows][c, :] .≤ (1 .- m[:c_w][c, :]) .* big_M)

    B = [br.b for br in r.branches]

    @constraint(m, [c in cases(r), edge in edge_ids(r)],
        m[:flows][c, edge] - sum(B[edge] * r.A'[bus, edge] * m[:ϕ][c, bus] for bus in buses(r)) ≤ big_M * m[:c_w][c, edge])
    @constraint(m, [c in cases(r), edge in edge_ids(r)],
        m[:flows][c, edge] - sum(B[edge] * r.A'[bus, edge] * m[:ϕ][c, bus] for bus in buses(r)) ≥ -big_M * m[:c_w][c, edge])
end

function OTS_flows_pst!(m::Model, r::TNR)
    big_M = 5 # TODO: change
    @info "HARDCODED OTS_flows_pst - bigM to $big_M"

    @constraint(m, [c in cases(r)], m[:γ_branch][c, :] .≤ m[:c_w][c, :] .* big_M)
    @constraint(m, [c in cases(r)], -m[:γ_branch][c, :] .≤ m[:c_w][c, :] .* big_M)

    @constraint(m, [c in cases(r)], m[:flows][c, :] .≤ (1 .- m[:c_w][c, :]) .* big_M)
    @constraint(m, [c in cases(r)], -m[:flows][c, :] .≤ (1 .- m[:c_w][c, :]) .* big_M)

    @constraint(m, [c in n_1cases(r)], m[:flows][c, r.outages[c]] == 0)  # nb_cases is the index of the base case

    D = spdiagm([e_index_for(r.g, e).b for e in edges(r.g)])
    B = r.A' * D * r.A
    N = nv(r.g)
    B_inv = inv(Matrix(B + fill(1 / N, N, N))) - fill(1 / N, N, N)
    p_to_f = D * r.A * B_inv
    δ_to_f = D - p_to_f * r.A' * D

    for c in cases(r)
        outage_phase = @expression(m, c == nb_cases(r) ? 0 : m[:β][c] .* [is_outage(r, c, e) for e in edge_ids(r)])
        @constraint(m, m[:flows][c, :] .== p_to_f * (m[:load][c, :] .- m[:gen][c, :]) .+ δ_to_f * (outage_phase .+ m[:γ_branch][c, :]))
    end
end

function TNR_flows_to_extra_bus!(m::Model, r::TNR)
    big_M = 5
    @info "HARDCODED TNR_flows_to_extra_bus - bigM to $big_M"

    for (bc_id, bus, sb) in ((bc_id, bc.bus, sb) for (bc_id, bc) in enumerate(r.bus_confs) for sb in bc.subBuses)
        @constraint(m, [c in cases(r)],
            sum(m[:flows][c, br] * r.A[br, bus] for br in sb.branch_ids) - (m[:hatload][c, bc_id] - m[:hatgen][c, bc_id]) ≤ (1 - m[:v_bus][bc_id]) * big_M)
        @constraint(m, [c in cases(r)],
            -(sum(m[:flows][c, br] * r.A[br, bus] for br in sb.branch_ids) - (m[:hatload][c, bc_id] - m[:hatgen][c, bc_id])) ≤ (1 - m[:v_bus][bc_id]) * big_M)
    end
end

function TNR_flows_phases!(m::Model, r::TNR)
    big_M = 100 # TODO: change
    @info "HARDCODED TNR_flows - bigM to $big_M"

    @constraint(m, [c in cases(r)], m[:ϕ][c, r.bus_orig_id] == 0)

    # flows in the branches
    @constraint(m, [c in cases(r)], m[:flows][c, :] .≤ (1 .- m[:c_w][c, :]) .* big_M)
    @constraint(m, [c in cases(r)], -m[:flows][c, :] .≤ (1 .- m[:c_w][c, :]) .* big_M)
    @constraint(m, [c in cases(r), edge in edge_ids(r)],
        m[:flows][c, edge] - r.branches[edge].b * (m[:ϕ_e][c, edge, 2] - m[:ϕ_e][c, edge, 1]) ≤ big_M * m[:c_w][c, edge])
    @constraint(m, [c in cases(r), edge in edge_ids(r)],
        m[:flows][c, edge] - r.branches[edge].b * (m[:ϕ_e][c, edge, 2] - m[:ϕ_e][c, edge, 1]) ≥ -big_M * m[:c_w][c, edge])

    # connect of the branch extremity phases to the bus phases
    for (bus, edge, c) in ((b,e,c) for b in buses(r), e in edge_ids(r), c in cases(r) if r.A[e,b] ≠ 0)
        ϕ_edge = @expression(m, r.A[edge, bus]==1 ? m[:ϕ_e][c, edge, 1] : m[:ϕ_e][c, edge, 2])
        if (bus, edge) in keys(r.bus_branch_to_conf_ids)
            bc = r.bus_branch_to_conf_ids[bus, edge][1] #TODO only one alternative conf per bus first !
            @constraint(m, m[:ϕ][c,bus] - ϕ_edge ≤  big_M * m[:v_bus][bc])
            @constraint(m, m[:ϕ][c,bus] - ϕ_edge ≥ -big_M * m[:v_bus][bc])
            @constraint(m, m[:ϕ_bc][c, bc] - ϕ_edge ≤  big_M * m[:v_bus][bc])
            @constraint(m, m[:ϕ_bc][c, bc] - ϕ_edge ≥ -big_M * m[:v_bus][bc])
        else
            @constraint(m, m[:ϕ][c,bus] - ϕ_edge == 0)
        end
    end

    # balance the flows in substations
    @constraint(m, [c in cases(r), bus in buses(r); bus ≠ r.bus_orig_id],
        m[:load][c, bus] - m[:gen][c, bus] == sum(r.A'[bus, e] * m[:flows][c, e] for e in edge_ids(r)))

        # balance flows in extra buses
    TNR_flows_to_extra_bus!(m,r)
end

function TNR_flows_pst!(m::Model, r::TNR)
    @info "TNR_flows_pst to be revisited with c_w rather than v_branch"
    big_M = 5 # TODO: change
    @info "HARDCODED TNR_flows - bigM to $big_M"

    if r.allow_branch_openings
        @constraint(m, [c in cases(r)], m[:γ_branch][c, :] .≤ m[:v_branch] .* big_M)
        @constraint(m, [c in cases(r)], -m[:γ_branch][c, :] .≤ m[:v_branch] .* big_M)

        @constraint(m, [c in cases(r)], m[:flows][c, :] .≤ (1 .- m[:v_branch]) .* big_M)
        @constraint(m, [c in cases(r)], -m[:flows][c, :] .≤ (1 .- m[:v_branch]) .* big_M)
    end

    @constraint(m, [c in cases(r), bc_id in bus_conf_ids(r)], m[:γ_bus][c, bc_id] ≤ m[:v_bus][bc_id] * big_M)
    @constraint(m, [c in cases(r), bc_id in bus_conf_ids(r)], -m[:γ_bus][c, bc_id] ≤ m[:v_bus][bc_id] * big_M)

    @constraint(m, [c in n_1cases(r)], m[:flows][c, r.outages[c]] == 0)  # nb_cases is the index of the base case

    D = spdiagm([e_index_for(r.g, e).b for e in edges(r.g)])
    B = r.A' * D * r.A
    N = nv(r.g)
    B_inv = inv(Matrix(B + fill(1 / N, N, N))) - fill(1 / N, N, N)
    p_to_f = D * r.A * B_inv
    δ_to_f = D - p_to_f * r.A' * D

    for c in cases(r)
        outage_phase = @expression(m, c == nb_cases(r) ? 0 : m[:β][c] .* [is_outage(r, c, e) for e in edge_ids(r)])
        line_opening_phase = @expression(m, r.allow_branch_openings ? m[:γ_branch][c, :] : zeros(edge_ids(r)))
        bus_splitting_phase = @expression(m, isempty(r.bus_confs) ? 0
                                             : [sum(br in sb.branch_ids ? (r.A[br, bc.bus] * m[:γ_bus][c, k]) : 0
                                                  for (k, bc) in enumerate(r.bus_confs) for sb in r.bus_confs[k].subBuses) for br = edge_ids(r)])

        @constraint(m, m[:flows][c, :] .== p_to_f * (m[:load][c, :] .- m[:gen][c, :])
                                           +
                                           δ_to_f * (outage_phase .+ line_opening_phase .+ bus_splitting_phase))
    end
    TNR_flows_to_extra_bus!(m,r)
end

function overload!(m::Model, r::TNR)
    p_max = [e_index_for(r.g, br).p_max for br in edges(r.g)]
    @constraint(m, [i in cases(r)], m[:flows][i, :] .≤ p_max)
    @constraint(m, [i in cases(r)], -m[:flows][i, :] .≤ p_max)
end

function c_w!(m::Model, r::TNR)
    if r.allow_branch_openings || r.OTS_only
        @constraint(m, [c in cases(r), e in edge_ids(r)], m[:c_w][c, e] ≥ is_outage(r, c, e))
        @constraint(m, [c in cases(r), e in edge_ids(r)], m[:c_w][c, e] ≥ m[:v_branch][e])
        @constraint(m, [c in cases(r), e in edge_ids(r)], m[:c_w][c, e] ≤ is_outage(r, c, e) + m[:v_branch][e])
    else
        @constraint(m, [c in cases(r), e in edge_ids(r)], m[:c_w][c, e] == is_outage(r, c, e))
    end
end

function OTS_N_connectedness!(m::Model, r::TNR)
    big_M2 = nb_buses(r)
    cases = (r.n_1_connectedness ? 1 : nb_cases(r)):nb_cases(r)        # nb_cases is the index of the base case

    bus = (e, s) -> findfirst(v -> r.A[e, v] == s, buses(r))

    # # c_flows are null when a line is open : c_w == 1
    @constraint(m, [c in cases, e in edge_ids(r)], m[:c_flows][c, e] ≤ big_M2 * (1 - m[:c_w][c, e]))
    @constraint(m, [c in cases, e in edge_ids(r)], -m[:c_flows][c, e] ≤ big_M2 * (1 - m[:c_w][c, e]))

    # # the first bus holds a generator with n_bus_tot - 1 - its number of possible confs,
    # # each other consumes 1 + (the number of possible confs on the bus)
    # # This ensures the c_flows are balanced at the substation level
    @constraint(m, [c in cases, bus in buses(r)],
        (bus == 1 ? -nb_buses(r) : 0) + 1 - sum(m[:c_flows][c, br] * r.A[br, bus] for br in edge_ids(r)) == 0)
end

function TNR_N_connectedness!(m::Model, r::TNR)
    nb_bus_tot = nb_buses(r) + (isempty(r.bus_confs) ? 0 : sum(length(r.bus_to_conf_ids[bus]) for bus in keys(r.bus_to_conf_ids)))
    big_M2 = nb_bus_tot + 1
    cases = (r.n_1_connectedness ? 1 : nb_cases(r)):nb_cases(r)        # nb_cases is the index of the base case

    # # c_flows are null when a line is open : c_w == 1
    @constraint(m, [c in cases, e in edge_ids(r)], m[:c_flows][c, e] ≤ big_M2 * (1 - m[:c_w][c, e]))
    @constraint(m, [c in cases, e in edge_ids(r)], -m[:c_flows][c, e] ≤ big_M2 * (1 - m[:c_w][c, e]))

    # # the first bus holds a generator with n_bus_tot - 1 - its number of possible confs,
    # # each other consumes 1 + (the number of possible confs on the bus)
    # # This ensures the c_flows are balanced at the substation level
    nb_conf = bus -> bus_nb_confs(r, bus)
    @constraint(m, [c in cases, bus in buses(r)],
        (bus == 1 ? -nb_bus_tot : 0) + 1 + nb_conf(bus) - sum(m[:c_flows][c, br] * r.A[br, bus] for br in edge_ids(r)) == 0)

    # But when confs are activated, the balance needs also to be ensured for the concerned subBuses
    for (bc_id, bc) in enumerate(r.bus_confs), c in cases
        c_flows_on_bus_bus = @expression(m,
            sum(m[:c_flows][c, br] * r.A[br, bc.bus] for sb in bc.subBuses for br in sb.branch_ids))
        @constraint(m, 1 - c_flows_on_bus_bus ≤ big_M2 * (1 - m[:v_bus][bc_id]))
        @constraint(m, -(1 - c_flows_on_bus_bus) ≤ big_M2 * (1 - m[:v_bus][bc_id]))
    end
end

function OTS_N_1_connectedness!(m::Model, r::TNR)
    @info "Bus orig is: $(r.bus_orig_id)"
    bigM_nb_v = nv(g) + 2

    @constraint(m, [c in n_1cases(r), e in edge_ids(r)], m[:cn1_flows][c, e] ≤ bigM_nb_v * (1 - m[:c_w][c, e]))
    @constraint(m, [c in n_1cases(r), e in edge_ids(r)], -m[:cn1_flows][c, e] ≤ bigM_nb_v * (1 - m[:c_w][c, e]))

    @constraint(m, [c in n_1cases(r), bus in [bus for bus in buses(r) if bus ≠ r.bus_orig_id]],
        m[:cn1_π][c, bus] - sum(r.A[e, bus] * m[:cn1_flows][c, e] for e in edge_ids(r)) == 0)
    @constraint(m, [c in n_1cases(r)],
        m[:cn1_p_orig][c] - sum(r.A[e, r.bus_orig_id] * m[:cn1_flows][c, e] for e in edge_ids(r)) == 0)

    @constraint(m, [c in n_1cases(r)], m[:cn1_π][c, r.bus_orig_id] == 1)

    @constraint(m, [c in n_1cases(r), bus in buses(r); bus ≠ r.bus_orig_id],
        m[:cn1_π][c, bus] ≤ sum(m[:cn1_ψ][c, bus, e] for e in incident(r, bus)))

    @constraint(m, [c in n_1cases(r), bus in buses(r), e in incident(r, bus)],
        m[:cn1_π][c, bus] ≥ m[:cn1_ψ][c, bus, e])

    @constraint(m, [c in n_1cases(r), bus in buses(r), e in incident(r, bus)],
        m[:cn1_ψ][c, bus, e] ≥ m[:cn1_π][c, opposite(r, e, bus)] - m[:c_w][c, e])
    @constraint(m, [c in n_1cases(r), bus in buses(r), e in incident(r, bus)],
        m[:cn1_ψ][c, bus, e] ≤ m[:cn1_π][c, opposite(r, e, bus)])
    @constraint(m, [c in n_1cases(r), bus in buses(r), e in incident(r, bus)],
        m[:cn1_ψ][c, bus, e] ≤ 1 - m[:c_w][c, e])
end

function TNR_N_1_connectednes!(m::Model, r::TNR)
    buses_wo_O = (k for k in buses(r) if k ≠ r.bus_orig_id)
    bigM_nb_v = nb_buses(r) + 2 + nb_bus_confs(r)

    for c in n_1cases(r), bus in buses(r), edge in incident(r, bus)
        @constraint(m, m[:cn1_χ][c, bus, edge] == m[:cn1_χ_circ][c, bus, edge]
                                                  +
                                                  sum(m[:cn1_hatπ][c, conf]
                                                      for conf in get(r.bus_to_conf_ids, opposite(r, edge, bus), Int[])
                                                      if edge in r.bus_confs[conf].subBuses[1].branch_ids))
        @constraint(m, m[:cn1_χ_circ][c, bus, edge] ≤ m[:cn1_π][c, opposite(r, edge, bus)])
        @constraint(m, m[:cn1_χ_circ][c, bus, edge] ≤ 1
                                                      -
                                                      sum(m[:v_bus][conf]
                                                          for conf in get(r.bus_to_conf_ids, opposite(r, edge, bus), Int[])
                                                          if edge in r.bus_confs[conf].subBuses[1].branch_ids))
        @constraint(m, m[:cn1_χ_circ][c, bus, edge] ≥ m[:cn1_π][c, opposite(r, edge, bus)]
                                                      -
                                                      sum(m[:v_bus][conf]
                                                          for conf in get(r.bus_to_conf_ids, opposite(r, edge, bus), Int[])
                                                          if edge in r.bus_confs[conf].subBuses[1].branch_ids))

        @constraint(m, m[:cn1_ψ][c, bus, edge] ≤ m[:cn1_χ][c, bus, edge])
        @constraint(m, m[:cn1_ψ][c, bus, edge] ≤ 1 - m[:c_w][c, edge])
        @constraint(m, m[:cn1_ψ][c, bus, edge] ≥ m[:cn1_χ][c, bus, edge] - m[:c_w][c, edge])
    end

    @constraint(m, [c in n_1cases(r),
            bus in buses(r),
            conf in get(r.bus_to_conf_ids, bus, Int[]),
            edge in r.bus_confs[conf].subBuses[1].branch_ids],
        m[:cn1_hat_𝚿][c, conf] ≥ m[:cn1_ψ][c, bus, edge])
    @constraint(m, [c in n_1cases(r),
            bus in buses(r),
            conf in get(r.bus_to_conf_ids, bus, Int[])],
        m[:cn1_hat_𝚿][c, conf] ≤
        sum(m[:cn1_ψ][c, bus, edge] for edge in r.bus_confs[conf].subBuses[1].branch_ids))

    @constraint(m, [c in n_1cases(r), conf in bus_conf_ids(r)],
        m[:cn1_hatπ][c, conf] ≤ m[:v_bus][conf])
    @constraint(m, [c in n_1cases(r), conf in bus_conf_ids(r)],
        m[:cn1_hatπ][c, conf] ≤ m[:cn1_hat_𝚿][c, conf])
    @constraint(m, [c in n_1cases(r), conf in bus_conf_ids(r)],
        m[:cn1_hatπ][c, conf] ≥ m[:v_bus][conf] + m[:cn1_hat_𝚿][c, conf] - 1)

    for ((bus, branch), confs) in r.bus_branch_to_conf_ids
        @constraint(m, [conf in confs], m[:cn1_U][bus, branch] ≥ m[:v_bus][conf])
        @constraint(m, m[:cn1_U][bus, branch] ≤ sum(m[:v_bus][conf] for conf in confs))
    end

    for bus in buses(r), edge in incident(r, bus)
        if (bus, edge) in keys(r.bus_branch_to_conf_ids)
            @constraint(m, [c in n_1cases(r)], m[:cn1_𝚿][c, bus, edge] ≤ 1 - m[:cn1_U][bus, edge])
            @constraint(m, [c in n_1cases(r)], m[:cn1_𝚿][c, bus, edge] ≤ m[:cn1_ψ][c, bus, edge])
            @constraint(m, [c in n_1cases(r)], m[:cn1_𝚿][c, bus, edge] ≥ m[:cn1_ψ][c, bus, edge] - m[:cn1_U][bus, edge])
        else
            @constraint(m, [c in n_1cases(r)], m[:cn1_𝚿][c, bus, edge] == m[:cn1_ψ][c, bus, edge])
        end
    end

    @constraint(m, [c in n_1cases(r), bus in buses(r), edge in incident(r, bus)],
        m[:cn1_π][c, bus] ≥ m[:cn1_𝚿][c, bus, edge])
    @constraint(m, [c in n_1cases(r), bus in buses(r)],
        m[:cn1_π][c, bus] ≤ sum(m[:cn1_𝚿][c, bus, edge] for edge in incident(r, bus)))

    @constraint(m, [c in n_1cases(r), e in edge_ids(r)], m[:cn1_flows][c, e] ≤ bigM_nb_v * (1 - m[:c_w][c, e]))
    @constraint(m, [c in n_1cases(r), e in edge_ids(r)], -m[:cn1_flows][c, e] ≤ bigM_nb_v * (1 - m[:c_w][c, e]))

    @constraint(m, [c in n_1cases(r), bus_conf_ids(r)], m[:cn1_hatπ][c, :] .≤ m[:v_bus][:])

    for (bc_id, bus, sb) in ((bc_id, bc.bus, sb) for (bc_id, bc) in enumerate(r.bus_confs) for sb in bc.subBuses)
        @constraint(m, [c in n_1cases(r)],
            sum(m[:cn1_flows][c, edge] * r.A[edge, bus] for edge in sb.branch_ids) - m[:cn1_hatπ][c, bc_id] ≤ (1 - m[:v_bus][bc_id]) * bigM_nb_v)
        @constraint(m, [c in n_1cases(r)],
            sum(m[:cn1_flows][c, edge] * r.A[edge, bus] for edge in sb.branch_ids) - m[:cn1_hatπ][c, bc_id] ≥ -(1 - m[:v_bus][bc_id]) * bigM_nb_v)
    end

    @constraint(m, [c in n_1cases(r), bus in buses_wo_O],
        sum(m[:cn1_hatπ][c, bc_id] for bc_id in get(r.bus_to_conf_ids, bus, Int[])) + m[:cn1_π][c, bus] - sum(r.A[edge, bus] * m[:cn1_flows][c, edge] for edge in incident(r, bus)) == 0)

    @constraint(m, [c in n_1cases(r)], m[:cn1_π][c, r.bus_orig_id] == 1)

    @constraint(m, [c in n_1cases(r)],
        m[:cn1_p_orig][c] + sum(r.A[edge, r.bus_orig_id] * m[:cn1_flows][c, edge] for edge in incident(r, r.bus_orig_id)) == 0)
end

function N_balance!(m::Model, r::TNR)
    for bus in buses(r)
        if r.p[bus] < 0
            @constraint(m, m[:load][nb_cases(r), bus] == 0)
            @constraint(m, m[:gen][nb_cases(r), bus] == -r.p[bus])
        else
            @constraint(m, m[:gen][nb_cases(r), bus] == 0)
            @constraint(m, m[:load][nb_cases(r), bus] == r.p[bus])
        end
    end
end

function N_balance_OPF!(m::Model, r::TNR)
    for bus in buses(r)
        if r.p[bus] < 0
            @constraint(m, m[:loadshed][bus] == 0)
            @constraint(m, m[:load][nb_cases(r), bus] == 0)
            @constraint(m, m[:gen][nb_cases(r), bus] == -r.p[bus])
        else
            @constraint(m, m[:gen][nb_cases(r), bus] == 0)
            @constraint(m, m[:loadshed][bus] ≤ r.p[bus])
            @constraint(m, m[:load][nb_cases(r), bus] == r.p[bus] - m[:loadshed][bus])
        end
    end
end

function OTS_balance!(m::Model, r::TNR)
    bigM = 10 # TODO: to define
    @info "HARDCODED OTS_balance - bigM to $bigM"

    m[:cn1_π] = m[:cn1_π]
    for c in n_1cases(r), bus in buses(r)
        if r.p[bus] < 0
            @constraint(m, m[:load][c, bus] == 0)
            @constraint(m, m[:gen][c, bus] + m[:σ][c] * r.p[bus] ≤ bigM * (1 - m[:cn1_π][c, bus]))
            @constraint(m, -(m[:gen][c, bus] + m[:σ][c] * r.p[bus]) ≤ bigM * (1 - m[:cn1_π][c, bus]))
            @constraint(m, m[:gen][c, bus] ≤ bigM * m[:cn1_π][c, bus])
            @constraint(m, -m[:gen][c, bus] ≤ bigM * m[:cn1_π][c, bus])
        else
            @constraint(m, m[:gen][c, bus] == 0)
            @constraint(m, m[:load][c, bus] - r.p[bus] ≤ bigM * (1 - m[:cn1_π][c, bus]))
            @constraint(m, -(m[:load][c, bus] - r.p[bus]) ≤ bigM * (1 - m[:cn1_π][c, bus]))
            @constraint(m, m[:load][c, bus] ≤ bigM * m[:cn1_π][c, bus])
            @constraint(m, -m[:load][c, bus] ≤ bigM * m[:cn1_π][c, bus])
        end
    end

    @constraint(m, [c in n_1cases(r)], sum(m[:load][c, :] .- m[:gen][c, :]) == 0)
end

function OTS_balance_OPF!(m::Model, r::TNR)
    bigM = 10 # TODO: to define
    @info "HARDCODED OTS_balance_OPF - bigM to $bigM"

    m[:cn1_π] = m[:cn1_π]
    for bus in buses(r)
        if r.p[bus] < 0
            @constraint(m, m[:loadshed][bus] == 0)
        else
            @constraint(m, m[:loadshed][bus] ≤ r.p[bus])
        end
        for c in cases(r)
            extended_π = @expression(m, c == n_case(r) ? 1 : m[:cn1_π][c, bus])
            if r.p[bus] < 0
                @constraint(m, m[:load][c, bus] == 0)
                @constraint(m, m[:gen][c, bus] + m[:σ][c] * r.p[bus] ≤ bigM * (1 - extended_π))
                @constraint(m, -(m[:gen][c, bus] + m[:σ][c] * r.p[bus]) ≤ bigM * (1 - extended_π))
                @constraint(m, m[:gen][c, bus] ≤ bigM * extended_π)
                @constraint(m, -m[:gen][c, bus] ≤ bigM * extended_π)
            else
                @constraint(m, m[:gen][c, bus] == 0)
                @constraint(m, m[:load][c, bus] + m[:loadshed][bus] - r.p[bus] ≤ bigM * (1 - extended_π))
                @constraint(m, -(m[:load][c, bus] + m[:loadshed][bus] - r.p[bus]) ≤ bigM * (1 - extended_π))
                @constraint(m, m[:load][c, bus] ≤ bigM * extended_π)
                @constraint(m, -m[:load][c, bus] ≤ bigM * extended_π)
            end
        end
    end

    @constraint(m, [c in cases(r)], sum(m[:load][c, :] .- m[:gen][c, :]) == 0)
end

function TNR_balance!(m::Model, r::TNR)
    bigM = 2 * maximum(abs, r.p)

    @variable(m, cn1_hatd_u[bus_conf_ids(r)] ≥ 0)
    @variable(m, cn1_hatg_u[n_1cases(r), bus_conf_ids(r)] ≥ 0)

    hatp = [bc.subBuses[1].p for bc in r.bus_confs]

    # N connectedness is supposed to be ensured. hat_p_s = u_s * p_conf_s
    for bc_id in bus_conf_ids(r)
        _p = hatp[bc_id]
        if _p > 0
            @constraint(m, m[:hatload][nb_cases(r), bc_id] ≤ bigM * m[:v_bus][bc_id])
            @constraint(m, m[:hatload][nb_cases(r), bc_id] - _p ≤ bigM * (1 - m[:v_bus][bc_id]))
            @constraint(m, -(m[:hatload][nb_cases(r), bc_id] - _p) ≤ bigM * (1 - m[:v_bus][bc_id]))
            @constraint(m, m[:hatgen][nb_cases(r), bc_id] == 0)
        elseif _p < 0
            @constraint(m, -m[:hatgen][nb_cases(r), bc_id] ≤ bigM * m[:v_bus][bc_id])
            @constraint(m, m[:hatgen][nb_cases(r), bc_id] + _p ≤ bigM * (1 - m[:v_bus][bc_id]))
            @constraint(m, -(m[:hatgen][nb_cases(r), bc_id] + _p) ≤ bigM * (1 - m[:v_bus][bc_id]))
            @constraint(m, m[:hatload][nb_cases(r), bc_id] == 0)
        else
            @constraint(m, m[:hatgen][nb_cases(r), bc_id] == 0)
            @constraint(m, m[:hatload][nb_cases(r), bc_id] == 0)
        end
    end

    for c in n_1cases(r)
        @constraint(m, sum(m[:gen][c, bus] - m[:load][c, bus] for bus in buses(r))
                       ==
                       0)
        for bc_id in bus_conf_ids(r)
            _p = hatp[bc_id]
            if _p > 0
                @constraint(m, m[:hatload][c, bc_id] ≤ bigM * m[:cn1_hatπ][c, bc_id])
                @constraint(m, m[:hatload][c, bc_id] - _p ≤ bigM * (1 - m[:cn1_hatπ][c, bc_id]))
                @constraint(m, -(m[:hatload][c, bc_id] - _p) ≤ bigM * (1 - m[:cn1_hatπ][c, bc_id]))
                @constraint(m, m[:hatgen][c, bc_id] == 0)
            elseif _p < 0
                # println("c:$c, bc_id:$bc_id hatgen:$hatgen")
                @constraint(m, m[:hatgen][c, bc_id] ≤ bigM * m[:cn1_hatπ][c, bc_id])
                @constraint(m, m[:hatgen][c, bc_id] + _p * m[:σ][c] ≤ bigM * (1 - m[:cn1_hatπ][c, bc_id]))
                @constraint(m, -(m[:hatgen][c, bc_id] + _p * m[:σ][c]) ≤ bigM * (1 - m[:cn1_hatπ][c, bc_id]))
                @constraint(m, m[:hatload][c, bc_id] == 0)
            else
                @constraint(m, m[:hatgen][c, bc_id] == 0)
                @constraint(m, m[:hatload][c, bc_id] == 0)
            end
        end
    end

    for bc_id in bus_conf_ids(r)
        _p = hatp[bc_id]
        if _p > 0
            @constraint(m, [c in n_1cases(r)], cn1_hatg_u[c, bc_id] == 0)
            @constraint(m, cn1_hatd_u[bc_id] ≤ bigM * m[:v_bus][bc_id])
            @constraint(m, cn1_hatd_u[bc_id] - _p ≤ bigM * (1 - m[:v_bus][bc_id]))
            @constraint(m, -(cn1_hatd_u[bc_id] - _p) ≤ bigM * (1 - m[:v_bus][bc_id]))
        elseif _p < 0
            @constraint(m, cn1_hatd_u[bc_id] == 0)
            for c in n_1cases(r)
                @constraint(m, cn1_hatg_u[c, bc_id] ≤ bigM * m[:v_bus][bc_id])
                @constraint(m, cn1_hatg_u[c, bc_id] + m[:σ][c] * _p ≤ bigM * (1 - m[:v_bus][bc_id]))
                @constraint(m, -(cn1_hatg_u[c, bc_id] + m[:σ][c] * _p) ≤ bigM * (1 - m[:v_bus][bc_id]))
            end
        else
            @constraint(m, [c in n_1cases(r)], cn1_hatg_u[c, bc_id] == 0)
            @constraint(m, cn1_hatd_u[bc_id] == 0)
        end
    end

    for c in n_1cases(r), bus in buses(r)
        bc_ids = get(r.bus_to_conf_ids, bus, Int[])
        @constraint(m, m[:load][c, bus] - sum(m[:hatload][c, bc_id] for bc_id in bc_ids) ≤ bigM * m[:cn1_π][c, bus])
        @constraint(m, -(m[:load][c, bus] - sum(m[:hatload][c, bc_id] for bc_id in bc_ids)) ≤ bigM * m[:cn1_π][c, bus])
        @constraint(m, m[:gen][c, bus] - sum(m[:hatgen][c, bc_id] for bc_id in bc_ids) ≤ bigM * m[:cn1_π][c, bus])
        @constraint(m, -(m[:gen][c, bus] - sum(m[:hatgen][c, bc_id] for bc_id in bc_ids)) ≤ bigM * m[:cn1_π][c, bus])
        @constraint(m, m[:load][c, bus] - maximum((0, r.p[bus]))
                       +
                       sum(cn1_hatd_u[bc_id] - m[:hatload][c, bc_id] for bc_id in bc_ids) ≤ bigM * (1 - m[:cn1_π][c, bus]))
        @constraint(m, -(m[:load][c, bus] - maximum((0, r.p[bus]))
                         +
                         sum(cn1_hatd_u[bc_id] - m[:hatload][c, bc_id] for bc_id in bc_ids)) ≤ bigM * (1 - m[:cn1_π][c, bus]))
        @constraint(m, m[:gen][c, bus] + m[:σ][c] * minimum((0, r.p[bus]))
                       + sum(cn1_hatg_u[c, bc_id] - m[:hatgen][c, bc_id] for bc_id in bc_ids) ≤ bigM * (1 - m[:cn1_π][c, bus]))
        @constraint(m, -(m[:gen][c, bus] + m[:σ][c] * minimum((0, r.p[bus]))
                         + sum(cn1_hatg_u[c, bc_id] - m[:hatgen][c, bc_id] for bc_id in bc_ids)) ≤ bigM * (1 - m[:cn1_π][c, bus]))
    end
end

function loadloss!(m::Model, r::TNR)
    @constraint(m, [c in n_1cases(r)], m[:lostload][c] == sum(r.p[bus] - m[:load][c, bus] for bus in buses(r) if r.p[bus] ≥ 0))
end

function secured_dc_OTS(g::MetaGraph;
    contingencies::AbstractArray{Int}=Int[],
    n_1_connectedness::Bool=false,
    bus_confs::Vector{BusConf}=BusConf[],
    allow_branch_openings::Bool=true,
    OTS_only::Bool=false,
    tnr_pf::TNR_PF_TYPE=tnr_pf_phase,
    opf::Bool,
    bus_orig::String)

    model, r = init_model(g, contingencies, n_1_connectedness, bus_confs, allow_branch_openings, OTS_only, tnr_pf, opf, bus_orig)
    create_variables!(model, r)

    c_w!(model, r)

    if !r.opf
        N_balance!(model, r)
    end

    if r.OTS_only
        if r.tnr_pf == tnr_pf_pst
            OTS_flows_pst!(model, r)
        elseif r.tnr_pf == tnr_pf_phase
            OTS_flows_phases!(model, r)
        end
        OTS_N_connectedness!(model, r)

        OTS_N_1_connectedness!(model, r)

        if r.opf
            OTS_balance_OPF!(model, r)
        else
            OTS_balance!(model, r)
        end

    else
        if r.tnr_pf == tnr_pf_pst
            TNR_flows_pst!(model, r)
        elseif r.tnr_pf == tnr_pf_phase
            TNR_flows_phases!(model, r)
        end
        TNR_N_connectedness!(model, r)
        TNR_N_1_connectednes!(model, r)
        TNR_balance!(model, r)

        @constraint(model, [bus in keys(r.bus_to_conf_ids)], sum(model[:v_bus][bc] for bc in r.bus_to_conf_ids[bus]) ≤ 1)
    end

    loadloss!(model, r)

    overload!(model, r)

    @objective(model, Min,
        sum(model[:lostload][c] for c in n_1cases(r)) +
        00 * (r.opf ? sum(model[:loadshed][bus] for bus in buses(r)) : 0))
    # + 0.01  * (allow_branch_openings ? sum(model[:v_branch]) : 0)) 
    # + 0.01 * (!isempty(bus_confs)   ? sum(model[:v_bus])    : 0) )
    model, r
end