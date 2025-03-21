using Logging
using JuMP
include("../src/tnr.jl")
include("equivalent.jl")

function init_model(r::TNR)
    # model = direct_model(Gurobi.Optimizer())
    model = Model(Gurobi.Optimizer)
    # set_silent(model)
    set_optimizer_attribute(model, "DualReductions", 0)
    log_path = joinpath("tmp", "my_log_file.txt")
    open(log_path, "a") do f
        write(f, "________________________________________________________________\n")
        write(f, string(r), "\n")
    end
    set_optimizer_attribute(model, "LogFile", log_path)
    model.ext[:r] = r
    model, r
end

function init_model(g::MetaGraph,
    contingencies::AbstractArray{Int},
    n_1_connectedness::Bool,
    bus_confs::Vector{BusConf},
    allow_branch_openings::Bool,
    OTS_only::Bool,
    tnr_pf::TNR_PF_TYPE,
    opf::Bool,
    bus_orig::String)

    r = TNR(g, contingencies, n_1_connectedness, bus_confs, allow_branch_openings, OTS_only, tnr_pf, opf, bus_orig)

    init_model(r)
end

function create_variables!(m::Model, r::TNR, fixed_buses::Vector{Int}=Int[], eq::Union{Equivalent,Nothing}=nothing)

    @variable(m, c_w[cases(r), edge_ids(r)], Bin)
    @variable(m, flows[cases(r), edge_ids(r)])
    @variable(m, overload)

    @variable(m, gen[cases(r), bus in buses(r); !(bus in fixed_buses)] ≥ 0)
    @variable(m, load[cases(r), bus in buses(r); !(bus in fixed_buses)] ≥ 0)
    @variable(m, lostload[n_1cases(r)])
    
    @variable(m, c_slacks[bus in (isempty(fixed_buses) ? Int[r.bus_orig_id] : fixed_buses)])
    @variable(m, c_flows[(r.n_1_connectedness ? 1 : nb_cases(r)):nb_cases(r), edge_ids(r)])

    @variable(m, cn1_flows[n_1cases(r), edge_ids(r)])
    @variable(m, cn1_π[n_1cases(r), buses(r)], Bin)
    @variable(m, cn1_ψ[n_1cases(r), bus in buses(r), incident(r, bus)], Bin)

    if !isnothing(eq)
        @variable(m, flows_e[cases(r), eq.buses])
        if eq.allow_connector_opening
            @variable(m, ϕ_e[cases(r), eq.buses])
            @variable(m, v_branch_e[eq.buses], Bin)
            @variable(m, cn1_ψ_e_bus[n_1cases(r), e_bus in eq.buses], Bin)
            @variable(m, cn1_ψ_e_cc[n_1cases(r), e_cc in 1:length(eq.cc), eq.cc[e_cc].buses], Bin)
        end
        @variable(m, p_e[cases(r), eq.buses])
        @variable(m, l_cc[n_1cases(r), 1:length(eq.cc)])


        @variable(m, c_flows_e[r.n_1_connectedness ? cases(r) : (n_case(r):n_case(r)), eq.buses])

        @variable(m, cn1_flows_e[n_1cases(r), eq.buses])
        @variable(m, cn1_π_e_cc[n_1cases(r), 1:length(eq.cc)], Bin)
    end

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

function OTS_flows_phases!(m::Model, r::TNR, fixed_buses::Vector{Int}, fixed_ϕ::Vector{Float64}, eq::Union{Equivalent,Nothing})
    big_M = 100 #TODO: to change  
    @info "HARDCODED OTS_flows_phases bigM to $big_M"

    @constraint(m, [c in cases(r), (i, f_bus) in enumerate(fixed_buses)], m[:ϕ][c, f_bus] == fixed_ϕ[i])

    @constraint(m, [c in cases(r)], m[:flows][c, :] .≤ (1 .- m[:c_w][c, :]) .* big_M)
    @constraint(m, [c in cases(r)], -m[:flows][c, :] .≤ (1 .- m[:c_w][c, :]) .* big_M)

    B = [br.b for br in r.branches]

    @constraint(m, [c in cases(r), edge in edge_ids(r)],
        m[:flows][c, edge] - sum(B[edge] * r.A'[bus, edge] * m[:ϕ][c, bus] for bus in buses(r)) ≤ big_M * m[:c_w][c, edge])
    @constraint(m, [c in cases(r), edge in edge_ids(r)],
        m[:flows][c, edge] - sum(B[edge] * r.A'[bus, edge] * m[:ϕ][c, bus] for bus in buses(r)) ≥ -big_M * m[:c_w][c, edge])

    if isnothing(eq)
        @constraint(m, [c in cases(r), bus in buses(r); !(bus in fixed_buses)],
            m[:load][c, bus] - m[:gen][c, bus] ==
            sum(r.A'[bus, e] * m[:flows][c, e] for e in edge_ids(r)))
    else
        @constraint(m, [c in cases(r), bus in buses(r); !(bus in fixed_buses)],
            m[:load][c, bus] - m[:gen][c, bus] ==
            sum(r.A'[bus, e] * m[:flows][c, e] for e in edge_ids(r)) -
            (bus in eq.buses ? m[:flows_e][c, bus] : 0))

        if eq.allow_connector_opening
            @constraint(m, [c in cases(r), bus in eq.buses], m[:ϕ][c, bus] - m[:ϕ_e][c, bus] .≤ big_M .* m[:v_branch_e][bus])
            @constraint(m, [c in cases(r), bus in eq.buses], m[:ϕ][c, bus] - m[:ϕ_e][c, bus] .≥ -big_M .* m[:v_branch_e][bus])
            @constraint(m, [c in cases(r), bus in eq.buses], m[:flows_e][c, bus] ≤ (1 - m[:v_branch_e][bus]) * big_M)
            @constraint(m, [c in cases(r), bus in eq.buses], m[:flows_e][c, bus] ≥ -(1 - m[:v_branch_e][bus]) * big_M)
            for (bus_id, bus) in enumerate(eq.buses)
                @constraint(m, [c in cases(r)],
                    (m[:flows_e][c, bus] - m[:p_e][c, bus] -
                     sum(eq.B[bus_id, bus_ϕ_id] * m[:ϕ_e][c, bus_ϕ] for (bus_ϕ_id, bus_ϕ) in enumerate(eq.buses))) ≤
                    m[:v_branch_e][bus] * big_M)
                @constraint(m, [c in cases(r)],
                    (m[:flows_e][c, bus] - m[:p_e][c, bus] -
                     sum(eq.B[bus_id, bus_ϕ_id] * m[:ϕ_e][c, bus_ϕ] for (bus_ϕ_id, bus_ϕ) in enumerate(eq.buses))) ≥
                    -m[:v_branch_e][bus] * big_M)
            end
        else
            for (bus_id, bus) in enumerate(eq.buses)
                @constraint(m, [c in cases(r)],
                    (m[:flows_e][c, bus] - m[:p_e][c, bus] ==
                     sum(eq.B[bus_id, bus_ϕ_id] * m[:ϕ][c, bus_ϕ] for (bus_ϕ_id, bus_ϕ) in enumerate(eq.buses))))
            end
        end
        @constraint(m, [c in cases(r), cc in eq.cc],
            sum(m[:flows_e][c, bus] for bus in cc.buses) ==
            sum(m[:p_e][c, bus] for bus in cc.buses))
    end
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
    for (bus, edge, c) in ((b, e, c) for b in buses(r), e in edge_ids(r), c in cases(r) if r.A[e, b] ≠ 0)
        ϕ_edge = @expression(m, r.A[edge, bus] == 1 ? m[:ϕ_e][c, edge, 1] : m[:ϕ_e][c, edge, 2])
        if (bus, edge) in keys(r.bus_branch_to_conf_ids)
            bc = r.bus_branch_to_conf_ids[bus, edge][1] #TODO only one alternative conf per bus first !
            @constraint(m, m[:ϕ][c, bus] - ϕ_edge ≤ big_M * m[:v_bus][bc])
            @constraint(m, m[:ϕ][c, bus] - ϕ_edge ≥ -big_M * m[:v_bus][bc])
            @constraint(m, m[:ϕ_bc][c, bc] - ϕ_edge ≤ big_M * m[:v_bus][bc])
            @constraint(m, m[:ϕ_bc][c, bc] - ϕ_edge ≥ -big_M * m[:v_bus][bc])
        else
            @constraint(m, m[:ϕ][c, bus] - ϕ_edge == 0)
        end
    end

    # balance the flows in substations
    @constraint(m, [c in cases(r), bus in buses(r); bus ≠ r.bus_orig_id],
        m[:load][c, bus] - m[:gen][c, bus] == sum(r.A'[bus, e] * m[:flows][c, e] for e in edge_ids(r)))

    # balance flows in extra buses
    TNR_flows_to_extra_bus!(m, r)
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
    TNR_flows_to_extra_bus!(m, r)
end

function overload!(m::Model, r::TNR, eq::Union{Equivalent,Nothing})
    p_max = [e_index_for(r.g, br).p_max for br in edges(r.g)]
    @constraint(m, [c in cases(r)], m[:flows][c, :] .≤ p_max)
    @constraint(m, [c in cases(r)], -m[:flows][c, :] .≤ p_max)

    if !isnothing(eq)
        if eq.allow_connector_opening
            @constraint(m, [c in cases(r), e_id in 1:length(eq.f)], eq.f[e_id] + sum(eq.F[e_id, :] .* m[:ϕ_e][c, :]) ≤ eq.f_max[e_id])
            @constraint(m, [c in cases(r), e_id in 1:length(eq.f)], eq.f[e_id] + sum(eq.F[e_id, :] .* m[:ϕ_e][c, :]) ≥ -eq.f_max[e_id])
        else
            @constraint(m, [c in cases(r), e_id in 1:length(eq.f)],
                eq.f[e_id] + sum(eq.F[e_id, bus_id] * m[:ϕ][c, bus] for (bus_id, bus) in enumerate(eq.buses)) ≤ eq.f_max[e_id])
            @constraint(m, [c in cases(r), e_id in 1:length(eq.f), bus in eq.buses],
                eq.f[e_id] + sum(eq.F[e_id, bus_id] .* m[:ϕ][c, bus] for (bus_id, bus) in enumerate(eq.buses)) ≥ -eq.f_max[e_id])
        end
    end
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

function OTS_N_connectedness!(m::Model, r::TNR, fixed_buses::Vector{Int}, eq::Union{Nothing,Equivalent})
    _fixed_buses = isempty(fixed_buses) ? Int[r.bus_orig_id] : fixed_buses
    big_M = nb_buses(r) + (isnothing(eq) ? 0 : length(eq.buses))
    cases = r.n_1_connectedness ? cases(r) : (n_case(r):n_case(r))

    @constraint(m, [c in cases, e in edge_ids(r)], m[:c_flows][c, e] ≤ big_M * (1 - m[:c_w][c, e]))
    @constraint(m, [c in cases, e in edge_ids(r)], -m[:c_flows][c, e] ≤ big_M * (1 - m[:c_w][c, e]))

    if isnothing(eq)
        @constraint(m, [c in cases, bus in buses(r); !(bus in _fixed_buses)],
            sum(m[:c_flows][c, br] * r.A[br, bus] for br in edge_ids(r)) == 1)
        @constraint(m, [c in cases, bus in _fixed_buses],
            sum(m[:c_flows][c, br] * r.A[br, bus] for br in edge_ids(r)) == m[:c_slacks][bus])
        @constraint(m, sum(m[:c_slacks][bus] for bus in _fixed_buses) == length(_fixed_buses) - nb_buses(r))

    else
        @constraint(m, [c in cases, bus in buses(r); !(bus in _fixed_buses)],
            sum(m[:c_flows][c, br] * r.A[br, bus] for br in edge_ids(r)) -
            ((bus in eq.buses) ? m[:c_flows_e][c, bus] : 0) == 1)
        @constraint(m, [c in cases, bus in _fixed_buses],
            sum(m[:c_flows][c, br] * r.A[br, bus] for br in edge_ids(r)) -
            ((bus in eq.buses) ? m[:c_flows_e][c, bus] : 0) == m[:c_slacks][bus])
        @constraint(m, sum(m[:c_slacks][bus] for bus in _fixed_buses) == length(_fixed_buses) - nb_buses(r) - length(eq.cc))

        if eq.allow_connector_opening
            @constraint(m, [c in cases, e_bus in eq.buses], m[:c_flows_e][c, e_bus] ≤ big_M * (1 - m[:v_branch_e][e_bus]))
            @constraint(m, [c in cases, e_bus in eq.buses], -m[:c_flows_e][c, e_bus] ≤ big_M * (1 - m[:v_branch_e][e_bus]))
        end

        @constraint(m, [c in cases, cc in eq.cc], sum(m[:c_flows_e][c, br] for br in cc.buses) == 1)
    end

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

function OTS_N_1_connectedness!(m::Model, r::TNR, fixed_buses::Vector{Int}, eq::Union{Nothing,Equivalent})
    bigM_nb_v = nb_buses(r) + 2 + (isnothing(eq) ? 0 : length(eq.buses))

    _fixed_buses = isempty(fixed_buses) ? Int[r.bus_orig_id] : fixed_buses

    # Set of constraints to ensure cn1_π is 0 in unenergized ares.
    @constraint(m, [c in n_1cases(r), e in edge_ids(r)], m[:cn1_flows][c, e] ≤ bigM_nb_v * (1 - m[:c_w][c, e]))
    @constraint(m, [c in n_1cases(r), e in edge_ids(r)], -m[:cn1_flows][c, e] ≤ bigM_nb_v * (1 - m[:c_w][c, e]))

    if isnothing(eq)
        @constraint(m, [c in n_1cases(r), bus in [bus for bus in buses(r) if !(bus in _fixed_buses)]],
            sum(r.A[e, bus] * m[:cn1_flows][c, e] for e in edge_ids(r)) ==
            m[:cn1_π][c, bus])
    else
        @constraint(m, [c in n_1cases(r), bus in [bus for bus in buses(r) if !(bus in _fixed_buses)]],
            sum(r.A[e, bus] * m[:cn1_flows][c, e] for e in edge_ids(r)) -
            (bus in eq.buses ? m[:cn1_flows_e][c, bus] : 0) ==
            m[:cn1_π][c, bus])
        if eq.allow_connector_opening
            @constraint(m, [c in n_1cases(r), bus in eq.buses], m[:cn1_flows_e][c, bus] ≤ bigM_nb_v * (1 - m[:v_branch_e][bus]))
            @constraint(m, [c in n_1cases(r), bus in eq.buses], -m[:cn1_flows_e][c, bus] ≤ bigM_nb_v * (1 - m[:v_branch_e][bus]))
        end
        @constraint(m, [c in n_1cases(r), cc_id in 1:length(eq.cc)],
            sum(m[:cn1_flows_e][c, e] for e in eq.cc[cc_id].buses) == m[:cn1_π_e_cc][c, cc_id])
    end


    # Set of constraints to ensure cn1_π is 1 in energized ares.
    @constraint(m, [c in n_1cases(r), fixed_bus in _fixed_buses], m[:cn1_π][c, fixed_bus] == 1)

    @constraint(m, [c in n_1cases(r), bus in buses(r), e in incident(r, bus); !(bus in _fixed_buses)],
        m[:cn1_π][c, bus] ≥ m[:cn1_ψ][c, bus, e])

    @constraint(m, [c in n_1cases(r), bus in buses(r), e in incident(r, bus); !(bus in _fixed_buses)],
        m[:cn1_ψ][c, bus, e] ≥ m[:cn1_π][c, opposite(r, e, bus)] - m[:c_w][c, e])
    @constraint(m, [c in n_1cases(r), bus in buses(r), e in incident(r, bus); !(bus in _fixed_buses)],
        m[:cn1_ψ][c, bus, e] ≤ m[:cn1_π][c, opposite(r, e, bus)])
    @constraint(m, [c in n_1cases(r), bus in buses(r), e in incident(r, bus); !(bus in _fixed_buses)],
        m[:cn1_ψ][c, bus, e] ≤ 1 - m[:c_w][c, e])

    if isnothing(eq)
        @constraint(m, [c in n_1cases(r), bus in buses(r); !(bus in _fixed_buses)],
            m[:cn1_π][c, bus] ≤ sum(m[:cn1_ψ][c, bus, e] for e in incident(r, bus)))
    else
        if eq.allow_connector_opening
            @constraint(m, [c in n_1cases(r), bus in buses(r); !(bus in _fixed_buses)],
                m[:cn1_π][c, bus] ≤
                sum(m[:cn1_ψ][c, bus, e] for e in incident(r, bus)) +
                (bus in eq.buses ? m[:cn1_ψ_e_bus][c, bus] : 0))
            @constraint(m, [c in n_1cases(r), bus in eq.buses],
                m[:cn1_π][c, bus] ≥ m[:cn1_ψ_e_bus][c, bus])

            @constraint(m, [c in n_1cases(r), bus_id in 1:length(eq.buses)],
                m[:cn1_ψ_e_bus][c, eq.buses[bus_id]] ≤ m[:cn1_π_e_cc][c, eq.cc_id[bus_id]])
            @constraint(m, [c in n_1cases(r), bus in eq.buses],
                m[:cn1_ψ_e_bus][c, bus] ≤ 1 - m[:v_branch_e][bus])
            @constraint(m, [c in n_1cases(r), bus_id in 1:length(eq.buses)],
                m[:cn1_ψ_e_bus][c, eq.buses[bus_id]] ≥ m[:cn1_π_e_cc][c, eq.cc_id[bus_id]] - m[:v_branch_e][eq.buses[bus_id]])

            @constraint(m, [c in n_1cases(r), cc_id in 1:length(eq.cc), e_bus in eq.cc[cc_id].buses],
                m[:cn1_π_e_cc][c, cc_id] ≥ m[:cn1_ψ_e_cc][c, cc_id, e_bus])
            @constraint(m, [c in n_1cases(r), cc_id in 1:length(eq.cc)],
                m[:cn1_π_e_cc][c, cc_id] ≤ sum(m[:cn1_ψ_e_cc][c, cc_id, e_bus] for e_bus in eq.cc[cc_id].buses))

            @constraint(m, [c in n_1cases(r), cc_id in 1:length(eq.cc), e_bus in eq.cc[cc_id].buses],
                m[:cn1_ψ_e_cc][c, cc_id, e_bus] ≤ m[:cn1_π][c, e_bus])
            @constraint(m, [c in n_1cases(r), cc_id in 1:length(eq.cc), e_bus in eq.cc[cc_id].buses],
                m[:cn1_ψ_e_cc][c, cc_id, e_bus] ≤ 1 - m[:v_branch_e][e_bus])
            @constraint(m, [c in n_1cases(r), cc_id in 1:length(eq.cc), e_bus in eq.cc[cc_id].buses],
                m[:cn1_ψ_e_cc][c, cc_id, e_bus] ≥ m[:cn1_π][c, e_bus] - m[:v_branch_e][e_bus])
        else
            @constraint(m, [c in n_1cases(r), bus in buses(r); !(bus in _fixed_buses)],
                m[:cn1_π][c, bus] ≤
                sum(m[:cn1_ψ][c, bus, e] for e in incident(r, bus)) +
                (bus in eq.buses ? m[:cn1_π_e_cc][c, get_cci(eq, bus)] : 0))
            @constraint(m, [c in n_1cases(r), bus in eq.buses],
                m[:cn1_π][c, bus] ≥ m[:cn1_π_e_cc][c, get_cci(eq, bus)])
        end
    end
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

end

function N_balance!(m::Model, r::TNR, fixed_buses::Vector{Int}, eq::Union{Nothing,Equivalent})
    for bus in buses(r)
        (bus in fixed_buses) && continue
        @constraint(m, m[:load][n_case(r), bus] == max(0, r.p[bus]))
        @constraint(m, m[:gen][n_case(r), bus] == -min(0, r.p[bus]))
    end
    if !isnothing(eq)
        for (i, bus) in enumerate(eq.buses)
            @constraint(m, m[:p_e][n_case(r), bus] == eq.p[i])
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

function OTS_balance!(m::Model, r::TNR, fixed_buses::Vector{Int}, eq::Union{Nothing,Equivalent})
    bigM = 10 # TODO: to define
    @info "HARDCODED OTS_balance - bigM to $bigM"

    for c in n_1cases(r), bus in buses(r)
        (bus in fixed_buses) && continue
        p = r.p[bus]
        if p < 0
            @constraint(m, m[:load][c, bus] == 0)
            @constraint(m, m[:gen][c, bus] + m[:σ][c] * p ≤ bigM * (1 - m[:cn1_π][c, bus]))
            @constraint(m, -(m[:gen][c, bus] + m[:σ][c] * p) ≤ bigM * (1 - m[:cn1_π][c, bus]))
            @constraint(m, m[:gen][c, bus] ≤ bigM * m[:cn1_π][c, bus])
            @constraint(m, -m[:gen][c, bus] ≤ bigM * m[:cn1_π][c, bus])
        else
            @constraint(m, m[:gen][c, bus] == 0)
            @constraint(m, m[:load][c, bus] - p ≤ bigM * (1 - m[:cn1_π][c, bus]))
            @constraint(m, -(m[:load][c, bus] - p) ≤ bigM * (1 - m[:cn1_π][c, bus]))
            @constraint(m, m[:load][c, bus] ≤ bigM * m[:cn1_π][c, bus])
            @constraint(m, -m[:load][c, bus] ≤ bigM * m[:cn1_π][c, bus])
        end
    end
    if !isnothing(eq)
        for c in n_1cases(r)
            for (i, bus) in enumerate(eq.buses)
                @constraint(m, m[:p_e][c, bus] - eq.p[i] - (m[:σ][c] - 1) * eq.sloped_g[i] ≤ bigM * m[:cn1_π_e_cc][c, eq.cc_id[i]])
                @constraint(m, m[:p_e][c, bus] - eq.p[i] - (m[:σ][c] - 1) * eq.sloped_g[i] ≥ -bigM * m[:cn1_π_e_cc][c, eq.cc_id[i]])
            end
            for (i, cc) in enumerate(eq.cc)
                @constraint(m, m[:l_cc][c, i] == cc.l * m[:cn1_π_e_cc][c, i])
            end
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

function loadloss!(m::Model, r::TNR, fixed_buses::Vector{Int}, eq::Union{Nothing,Equivalent})
    if isnothing(eq)
        @constraint(m, [c in n_1cases(r)],
            m[:lostload][c] == sum(r.p[bus] - m[:load][c, bus] for bus in buses(r) if !(bus in fixed_buses) && r.p[bus] ≥ 0))
    else
        @constraint(m, [c in n_1cases(r)],
            m[:lostload][c] ==
            sum(r.p[bus] - m[:load][c, bus] for bus in buses(r) if !(bus in fixed_buses) && r.p[bus] ≥ 0) +
            sum(eq.cc[cc_id].l - m[:l_cc][c, cc_id] for cc_id in 1:length(eq.cc)))
        # sum(eq.p[i] - m[:p_e][c, bus] for (i, bus) in enumerate(eq.buses) if eq.p[i] ≥ 0))
    end
end

function branch_status!(m, r::TNR, branch_status::Union{AbstractArray{Bool},Dict{Int,Any},AbstractArray{Int},AbstractArray{Tuple{String,String}},Nothing})
    isnothing(branch_status) && return
    if isa(branch_status, AbstractArray{Bool})
        @constraint(m, m[:v_branch] .== branch_status)
    elseif isa(branch_status, Dict{Int,Any})
        @constraint(m, [br in keys(branch_status)], m[:v_branch][br] == branch_status[br])
    elseif isa(branch_status, AbstractArray{Int})
        @constraint(m, m[:v_branch] .== [e in final_outages for e in edge_ids(r)])
    elseif isa(branch_status, AbstractArray{Tuple{String,String}})
        @constraint(m, m[:v_branch] .== [e_label_for(r.g, e) in final_outages for e in edge_ids(r)])
    end
end

function secured_dc_OTS(g::MetaGraph;
    contingencies::AbstractArray{Int}=Int[],
    n_1_connectedness::Bool=false,
    bus_confs::Vector{BusConf}=BusConf[],
    allow_branch_openings::Bool=true,
    OTS_only::Bool=false,
    tnr_pf::TNR_PF_TYPE=tnr_pf_phase,
    opf::Bool=false,
    bus_orig::String,
    fixed_buses::Union{Vector{String},Vector{Int}}=Int[],
    fixed_phases::Vector{Float64}=Float64[],
    branch_status::Union{AbstractArray{Bool},Dict{Int,Any},AbstractArray{Int},AbstractArray{Tuple{String,String}},Nothing}=nothing,
    eq::Union{Equivalent,Nothing}=nothing)

    model, r = init_model(g, contingencies, n_1_connectedness, bus_confs, allow_branch_openings, OTS_only, tnr_pf, opf, bus_orig)

    _fixed_phases = fixed_phases
    _fixed_buses = isa(fixed_buses, Vector{Int}) ? fixed_buses : to_code(g, fixed_buses)

    @info "Fixed bus: $_fixed_buses"
    create_variables!(model, r, _fixed_buses, eq)

    c_w!(model, r)

    if !r.opf
        N_balance!(model, r, _fixed_buses, eq)
    end

    if r.OTS_only
        if r.tnr_pf == tnr_pf_pst
            OTS_flows_pst!(model, r)
        elseif r.tnr_pf == tnr_pf_phase
            OTS_flows_phases!(model, r, _fixed_buses, _fixed_phases, eq)
        end
        OTS_N_connectedness!(model, r, _fixed_buses, eq)

        OTS_N_1_connectedness!(model, r, _fixed_buses, eq)

        if r.opf
            OTS_balance_OPF!(model, r)
        else
            OTS_balance!(model, r, _fixed_buses, eq)
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

    loadloss!(model, r, _fixed_buses, eq)

    overload!(model, r, eq)

    branch_status!(model, r, branch_status)

    @objective(model, Min,
        sum(model[:lostload][c] for c in n_1cases(r)) +
        00 * (r.opf ? sum(model[:loadshed][bus] for bus in buses(r)) : 0))
    # + 0.01  * (allow_branch_openings ? sum(model[:v_branch]) : 0)) 
    # + 0.01 * (!isempty(bus_confs)   ? sum(model[:v_bus])    : 0) )
    optimize!(model)
    !is_solved_and_feasible(model) && @warn "NOT FEASIBLE"
    model, r
end
