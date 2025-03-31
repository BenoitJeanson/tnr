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
    contingencies::AbstractArray{ELabel},
    n_1_connectedness::Bool,
    bus_confs::Vector{BusConf},
    allow_branch_openings::Bool,
    OTS_only::Bool,
    tnr_pf::TNR_PF_TYPE,
    opf::Bool,
    bus_orig::VLabel)

    r = TNR(g, contingencies, n_1_connectedness, bus_confs, allow_branch_openings, OTS_only, tnr_pf, opf, bus_orig)

    init_model(r)
end

function create_variables!(m::Model, r::TNR, fixed_buses::Vector{VLabel}, eq::Vector{Equivalent})

    @variable(m, c_w[cases(r), busfrom in labels(r.g), outneighbor_labels(r.g, busfrom)], Bin)
    @variable(m, flows[cases(r), busfrom in labels(r.g), outneighbor_labels(r.g, busfrom)])
    @variable(m, overload)

    @variable(m, gen[cases(r), bus in buses(r); !(bus in fixed_buses)] ≥ 0)
    @variable(m, load[cases(r), bus in buses(r); !(bus in fixed_buses)] ≥ 0)
    @variable(m, lostload[n_1cases(r)])

    @variable(m, c_slacks[bus in (isempty(fixed_buses) ? VLabel[r.bus_orig] : fixed_buses)])
    @variable(m, c_flows[(r.n_1_connectedness ? cases(r) : [n_case(r)]), busfrom in labels(r.g), outneighbor_labels(r.g, busfrom)])

    @variable(m, cn1_flows[n_1cases(r), busfrom in labels(r.g), outneighbor_labels(r.g, busfrom)])
    @variable(m, cn1_π[n_1cases(r), buses(r)], Bin)
    @variable(m, cn1_ψ[n_1cases(r), bus in buses(r), busfrom in labels(r.g), outneighbor_labels(r.g, busfrom)], Bin)

    @variable(m, v_branch[busfrom in labels(r.g), outneighbor_labels(r.g, busfrom)], Bin)

    @variable(m, ϕ[cases(r), buses(r)])

    if r.opf
        @variable(m, loadshed[buses(r)] ≥ 0)
        @variable(m, σ[cases(r)] ≥ 0)
    else
        @variable(m, σ[n_1cases(r)] ≥ 0)
    end

    # Equivalent
    @variable(m, flows_e[cases(r), q in eachindex(eq), eq[q].buses])
    @variable(m, inner_flows_e[cases(r), q in eachindex(eq), [str(e) for e in eq[q].internal_branches]])
    @variable(m, p_e[cases(r), q in eachindex(eq), eq[q].buses])
    @variable(m, l_cc[n_1cases(r), q in eachindex(eq), eachindex(eq[q].cc)])

    @variable(m, c_flows_e[r.n_1_connectedness ? cases(r) : [n_case(r)], q in eachindex(eq), eq[q].buses])

    @variable(m, cn1_flows_e[n_1cases(r), q in eachindex(eq), eq[q].buses])
    @variable(m, cn1_π_e_cc[n_1cases(r), q in eachindex(eq), eachindex(eq[q].cc)], Bin)

    # allow_connector_opening
    @variable(m, ϕ_e[cases(r), q in eachindex(eq), eq[q].buses; eq[q].allow_connector_opening])
    @variable(m, v_branch_e[q in eachindex(eq), eq[q].buses; eq[q].allow_connector_opening], Bin)
    @variable(m, cn1_ψ_e_bus[n_1cases(r), q in eachindex(eq), e_bus in eq[q].buses; eq[q].allow_connector_opening], Bin)
    @variable(m, cn1_ψ_e_cc[n_1cases(r), q in eachindex(eq), e_cc in eachindex(eq[q].cc), eq[q].cc[e_cc].buses; eq[q].allow_connector_opening], Bin)
end

function OTS_flows_phases!(m::Model, r::TNR, fixed_buses::Vector{VLabel}, fixed_ϕ::Vector{Float64}, eq::Vector{Equivalent})
    big_M = 100 #TODO: to change  
    @info "HARDCODED OTS_flows_phases bigM to $big_M"

    for c in cases(r)
        @constraint(m, [(i, f_bus) in enumerate(fixed_buses)], m[:ϕ][c, f_bus] == fixed_ϕ[i])

        @constraint(m, m[:flows][c, :, :] .≤ (1 .- m[:c_w][c, :, :]) .* big_M)
        @constraint(m, -m[:flows][c, :, :] .≤ (1 .- m[:c_w][c, :, :]) .* big_M)

        @constraint(m, [b in branches(r)],
            m[:flows][c, b...] - sum(r.g[b...].b * (m[:ϕ][c, to(b)] - m[:ϕ][c, from(b)])) ≤ big_M * m[:c_w][c, b...])
        @constraint(m, [b in branches(r)],
            m[:flows][c, b...] - sum(r.g[b...].b * (m[:ϕ][c, to(b)] - m[:ϕ][c, from(b)])) ≥ -big_M * m[:c_w][c, b...])

        if isnothing(eq)
            @constraint(m, [bus in buses(r); !(bus in fixed_buses)],
                m[:load][c, bus] - m[:gen][c, bus] ==
                sum(sign * m[:flows][c, b...] for (b, sign) in incident_signed(r, bus)))
        else
            @constraint(m, [bus in buses(r); !(bus in fixed_buses)],
                m[:load][c, bus] - m[:gen][c, bus] ==
                sum(sign * m[:flows][c, b...] for (b, sign) in incident_signed(r, bus)) +
                sum(bus in eq[q].buses ? m[:flows_e][c, q, bus] : 0 for q in eachindex(eq)))

            for q in eachindex(eq)
                if eq[q].allow_connector_opening
                    @constraint(m, [bus in eq[q].buses], m[:ϕ][c, bus] - m[:ϕ_e][c, q, bus] .≤ big_M .* m[:v_branch_e][q, bus])
                    @constraint(m, [bus in eq.buses], m[:ϕ][c, bus] - m[:ϕ_e][c, q, bus] .≥ -big_M .* m[:v_branch_e][q, bus])
                    @constraint(m, [bus in eq.buses], m[:flows_e][c, q, bus] ≤ (1 - m[:v_branch_e][q, bus]) * big_M)
                    @constraint(m, [bus in eq.buses], m[:flows_e][c, q, bus] ≥ -(1 - m[:v_branch_e][q, bus]) * big_M)
                    for (bus_id, bus) in enumerate(eq[q].buses)
                        @constraint(m,
                            (m[:flows_e][c, q, bus] - m[:p_e][c, q, bus] -
                             sum(eq[q].B[bus_id, bus_ϕ_id] * m[:ϕ_e][c, q, bus_ϕ] for (bus_ϕ_id, bus_ϕ) in enumerate(eq[q].buses))) ≤
                            m[:v_branch_e][q, bus] * big_M)
                        @constraint(m,
                            (m[:flows_e][c, q, bus] - m[:p_e][c, q, bus] -
                             sum(eq[q].B[bus_id, bus_ϕ_id] * m[:ϕ_e][c, q, bus_ϕ] for (bus_ϕ_id, bus_ϕ) in enumerate(eq[q].buses))) ≥
                            -m[:v_branch_e][q, bus] * big_M)
                    end
                else
                    for (bus_id, bus) in enumerate(eq[q].buses)
                        @constraint(m,
                            m[:flows_e][c, q, bus] - m[:p_e][c, q, bus] ==
                            sum(eq[q].B[bus_id, bus_ϕ_id] * m[:ϕ][c, bus_ϕ] for (bus_ϕ_id, bus_ϕ) in enumerate(eq[q].buses)))
                    end
                end
                @constraint(m, [cc in eq[q].cc],
                    sum(m[:flows_e][c, q, bus] for bus in cc.buses) ==
                    sum(m[:p_e][c, q, bus] for bus in cc.buses))
            end
        end
    end
end

function overload!(m::Model, r::TNR, eq::Vector{Equivalent})
    for c in cases(r)
        @constraint(m, [b in branches(r)], m[:flows][c, b...] .≤ r.g[b...].p_max)
        @constraint(m, [b in branches(r)], -m[:flows][c, b...] .≤ r.g[b...].p_max)

        for q in eachindex(eq)
            if eq[q].allow_connector_opening
                @constraint(m, [q in eachindex(eq), e_id in eachindex(eq[q].f)],
                    eq[q].f[e_id] + sum(eq[q].F[e_id, :] .* m[:ϕ_e][c, q, :]) ≤ eq[q].f_max[e_id])
                @constraint(m, [q in eachindex(eq), e_id in eachindex(eq.f)],
                    eq[q].f[e_id] + sum(eq[q].F[e_id, :] .* m[:ϕ_e][c, q, :]) ≥ -eq[q].f_max[e_id])
            else
                for q in eachindex(eq)
                    ib = [str(eq[q].internal_branches[e_id]) for e_id in eachindex(eq[q].f)]
                    @constraint(m, [e_id in eachindex(eq[q].f)],
                        m[:inner_flows_e][c, q, ib[e_id]] == eq[q].f[e_id] + sum(eq[q].F[e_id, bus_id] * m[:ϕ][c, bus] for (bus_id, bus) in enumerate(eq[q].buses)))
                    @constraint(m, [e_id in eachindex(eq[q].f)],
                        m[:inner_flows_e][c, q, ib[e_id]] ≤ eq[q].f_max[e_id])
                    @constraint(m, [e_id in eachindex(eq[q].f), bus in eq[q].buses],
                        m[:inner_flows_e][c, q, ib[e_id]] ≥ -eq[q].f_max[e_id])
                end
            end
        end
    end
end

function c_w!(m::Model, r::TNR)
    if r.allow_branch_openings || r.OTS_only
        @constraint(m, [c in cases(r), e in branches(r)], m[:c_w][c, e...] ≥ is_outage(c, e))
        @constraint(m, [c in cases(r), e in branches(r)], m[:c_w][c, e...] ≥ m[:v_branch][e...])
        @constraint(m, [c in cases(r), e in branches(r)], m[:c_w][c, e...] ≤ is_outage(c, e) + m[:v_branch][e...])
    else
        @constraint(m, [c in cases(r), e in branches(r)], m[:c_w][c, e...] == is_outage(c, e))
    end
end

function OTS_N_connectedness!(m::Model, r::TNR, fixed_buses::Vector{VLabel}, eq::Vector{Equivalent})
    _fixed_buses = isempty(fixed_buses) ? VLabel[r.bus_orig] : fixed_buses
    big_M = nb_buses(r) + (isempty(eq) ? 0 : sum(length(q.buses) for q in eq))
    cases = r.n_1_connectedness ? cases(r) : [n_case(r)]

    for c in cases
        @constraint(m, m[:c_flows][c, :, :] .≤ big_M .* (1 .- m[:c_w][c, :, :]))
        @constraint(m, -m[:c_flows][c, :, :] .≤ big_M .* (1 .- m[:c_w][c, :, :]))

        if isempty(eq)
            @constraint(m, [bus in buses(r); !(bus in _fixed_buses)],
                sum(sign * m[:c_flows][c, b...] for (b, sign) in incident_signed(r, bus)) == 1)
            @constraint(m, [bus in _fixed_buses],
                sum(sign * m[:c_flows][c, b...] for (b, sign) in incident_signed(r, bus)) == m[:c_slacks][bus])
            @constraint(m, sum(m[:c_slacks][bus] for bus in _fixed_buses) == length(_fixed_buses) - nb_buses(r))

        else
            @constraint(m, [bus in buses(r); !(bus in _fixed_buses)],
                sum(sign * m[:c_flows][c, b...] for (b, sign) in incident_signed(r, bus)) -
                sum((bus in eq[q].buses) ? m[:c_flows_e][c, q, bus] : 0 for q in eachindex(eq)) == 1)
            @constraint(m, [bus in _fixed_buses],
                sum(sign * m[:c_flows][c, b...] for (b, sign) in incident_signed(r, bus)) -
                sum((bus in eq[q].buses) ? m[:c_flows_e][c, q, bus] : 0 for q in eachindex(eq)) == m[:c_slacks][bus])
            @constraint(m, sum(m[:c_slacks][bus] for bus in _fixed_buses) == length(_fixed_buses) - nb_buses(r) - sum(length(q.cc) for q in eq))

            for q in eachindex(eq)
                if eq[q].allow_connector_opening
                    @constraint(m, [e_bus in eq.buses], m[:c_flows_e][c, q, e_bus] ≤ big_M * (1 - m[:v_branch_e][q, e_bus]))
                    @constraint(m, [e_bus in eq.buses], -m[:c_flows_e][c, q, e_bus] ≤ big_M * (1 - m[:v_branch_e][q, e_bus]))
                end
                @constraint(m, [cc in eq[q].cc], sum(m[:c_flows_e][c, q, br] for br in cc.buses) == 1)
            end

        end
    end

end

function OTS_N_1_connectedness!(m::Model, r::TNR, fixed_buses::Vector{VLabel}, eq::Vector{Equivalent})
    bigM_nb_v = nb_buses(r) + 2 + (isempty(eq) ? 0 : sum(length(q.buses) for q in eq))

    _fixed_buses = isempty(fixed_buses) ? VLabel[r.bus_orig] : fixed_buses

    for c in n_1cases(r)
        # Set of constraints to ensure cn1_π is 0 in unenergized ares.
        @constraint(m, [b in branches(r)], m[:cn1_flows][c, b...] .≤ bigM_nb_v .* (1 .- m[:c_w][c, b...]))
        @constraint(m, [e in branches(r)], -m[:cn1_flows][c, :, :] .≤ bigM_nb_v .* (1 .- m[:c_w][c, :, :]))

        if isempty(eq)
            @constraint(m, [bus in [bus for bus in buses(r) if !(bus in _fixed_buses)]],
                sum(sign * m[:cn1_flows][c, b...] for (b, sign) in incident_signed(r, bus)) ==
                m[:cn1_π][c, bus])
        else
            @constraint(m, [bus in [bus for bus in buses(r) if !(bus in _fixed_buses)]],
                sum(sign * m[:cn1_flows][c, b...] for (b, sign) in incident_signed(r, bus)) -
                sum(bus in eq[q].buses ? m[:cn1_flows_e][c, q, bus] : 0 for q in eachindex(eq)) ==
                m[:cn1_π][c, bus])
            for q in eachindex(eq)
                if eq[q].allow_connector_opening
                    @constraint(m, [bus in eq.buses], m[:cn1_flows_e][c, q, bus] ≤ bigM_nb_v * (1 - m[:v_branch_e][q, bus]))
                    @constraint(m, [bus in eq.buses], -m[:cn1_flows_e][c, q, bus] ≤ bigM_nb_v * (1 - m[:v_branch_e][q, bus]))
                end
                @constraint(m, [cc_id in eachindex(eq[q].cc)],
                    sum(m[:cn1_flows_e][c, q, e] for e in eq[q].cc[cc_id].buses) == m[:cn1_π_e_cc][c, q, cc_id])
            end
        end


        # Set of constraints to ensure cn1_π is 1 in energized ares.
        @constraint(m, [fixed_bus in _fixed_buses], m[:cn1_π][c, fixed_bus] == 1)

        @constraint(m, [bus in buses(r), b in incident(r, bus); !(bus in _fixed_buses)],
            m[:cn1_π][c, bus] ≥ m[:cn1_ψ][c, bus, b...])

        @constraint(m, [bus in buses(r), b in incident(r, bus); !(bus in _fixed_buses)],
            m[:cn1_ψ][c, bus, b...] ≥ m[:cn1_π][c, opposite(b, bus)] - m[:c_w][c, b...])
        @constraint(m, [bus in buses(r), b in incident(r, bus); !(bus in _fixed_buses)],
            m[:cn1_ψ][c, bus, b...] ≤ m[:cn1_π][c, opposite(b, bus)])
        @constraint(m, [bus in buses(r), b in incident(r, bus); !(bus in _fixed_buses)],
            m[:cn1_ψ][c, bus, b...] ≤ 1 - m[:c_w][c, b...])

        if isempty(eq)
            @constraint(m, [bus in buses(r); !(bus in _fixed_buses)],
                m[:cn1_π][c, bus] ≤ sum(m[:cn1_ψ][c, bus, b...] for b in incident(r, bus)))
        else
            for q in eachindex(eq)
                if eq[q].allow_connector_opening
                    @constraint(m, [bus in buses(r); !(bus in _fixed_buses)],
                        m[:cn1_π][c, bus] ≤
                        sum(m[:cn1_ψ][c, bus, b...] for b in incident(r, bus)) +
                        (bus in eq[q].buses ? m[:cn1_ψ_e_bus][c, q, bus] : 0))
                    @constraint(m, [bus in eq[q].buses],
                        m[:cn1_π][c, bus] ≥ m[:cn1_ψ_e_bus][c, q, bus])

                    @constraint(m, [bus_id in eachindex(eq[q].buses)],
                        m[:cn1_ψ_e_bus][c, q, eq[q].buses[bus_id]] ≤ m[:cn1_π_e_cc][c, q, eq[q].cc_id[bus_id]])
                    @constraint(m, [bus in eq[q].buses],
                        m[:cn1_ψ_e_bus][c, q, bus] ≤ 1 - m[:v_branch_e][q, bus])
                    @constraint(m, [bus_id in eachindex(eq[q].buses)],
                        m[:cn1_ψ_e_bus][c, q, eq[q].buses[bus_id]] ≥ m[:cn1_π_e_cc][c, q, eq[q].cc_id[bus_id]] - m[:v_branch_e][q, eq[q].buses[bus_id]])

                    @constraint(m, [cc_id in eachindex(eq.cc), e_bus in eq.cc[cc_id].buses],
                        m[:cn1_π_e_cc][c, q, cc_id] ≥ m[:cn1_ψ_e_cc][c, q, cc_id, e_bus])
                    @constraint(m, [cc_id in eachindex(eq.cc)],
                        m[:cn1_π_e_cc][c, q, cc_id] ≤ sum(m[:cn1_ψ_e_cc][c, q, cc_id, e_bus] for e_bus in eq[q].cc[cc_id].buses))

                    @constraint(m, [cc_id in eachindex(eq[q].cc), e_bus in eq[q].cc[cc_id].buses],
                        m[:cn1_ψ_e_cc][c, q, cc_id, e_bus] ≤ m[:cn1_π][c, e_bus])
                    @constraint(m, [cc_id in eachindex(eq[q].cc), e_bus in eq[q].cc[cc_id].buses],
                        m[:cn1_ψ_e_cc][c, q, cc_id, e_bus] ≤ 1 - m[:v_branch_e][q, e_bus])
                    @constraint(m, [cc_id in eachindex(eq[q].cc), e_bus in eq[q].cc[cc_id].buses],
                        m[:cn1_ψ_e_cc][c, q, cc_id, e_bus] ≥ m[:cn1_π][c, e_bus] - m[:v_branch_e][q, e_bus])
                else
                    @constraint(m, [bus in buses(r); !(bus in _fixed_buses)],
                        m[:cn1_π][c, bus] ≤
                        sum(m[:cn1_ψ][c, bus, br...] for br in incident(r, bus)) +
                        (bus in eq[q].buses ? m[:cn1_π_e_cc][c, q, get_cci(eq[q], bus)] : 0))
                    @constraint(m, [bus in eq[q].buses],
                        m[:cn1_π][c, bus] == m[:cn1_π_e_cc][c, q, get_cci(eq[q], bus)])
                end
            end
        end
    end
end

function N_balance!(m::Model, r::TNR, fixed_buses::Vector{VLabel}, eq::Vector{Equivalent})
    for bus in buses(r)
        (bus in fixed_buses) && continue
        @constraint(m, m[:load][n_case(r), bus] == max(0, r.g[bus]))
        @constraint(m, m[:gen][n_case(r), bus] == -min(0, r.g[bus]))
    end
    if !isnothing(eq)
        for q in eachindex(eq), (i, bus) in enumerate(eq[q].buses)
            @constraint(m, m[:p_e][n_case(r), q, bus] == eq[q].p[i])
        end
    end
end

function N_balance_OPF!(m::Model, r::TNR)
    for bus in buses(r)
        if r.g[bus] < 0
            @constraint(m, m[:loadshed][bus] == 0)
            @constraint(m, m[:load][n_case(r), bus] == 0)
            @constraint(m, m[:gen][n_cases(r), bus] == -r.g[bus])
        else
            @constraint(m, m[:gen][n_cases(r), bus] == 0)
            @constraint(m, m[:loadshed][bus] ≤ r.g[bus])
            @constraint(m, m[:load][n_cases(r), bus] == r.g[bus] - m[:loadshed][bus])
        end
    end
end

function OTS_balance!(m::Model, r::TNR, fixed_buses::Vector{VLabel}, eq::Vector{Equivalent})
    bigM = 10 # TODO: to define
    @info "HARDCODED OTS_balance - bigM to $bigM"

    if isempty(fixed_buses)
        if isempty(eq)
            @constraint(m, [c in n_1cases(r)], sum(m[:load][c, :] .- m[:gen][c, :]) == 0)
        else
            @constraint(m, [c in n_1cases(r)], sum(m[:load][c, :] .- m[:gen][c, :]) - sum(m[:p_e][c, :, :]) == 0)
        end
    else
        @constraint(m, m[:σ][:] .== 1)
    end

    for c in n_1cases(r), bus in buses(r)
        (bus in fixed_buses) && continue
        p = r.g[bus]
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
        if !isempty(eq)
            for q in eachindex(eq), (i, bus) in enumerate(eq[q].buses)
                @constraint(m, m[:p_e][c, q, bus] - eq[q].p[i] - (m[:σ][c] - 1) * eq[q].gen_slope[i] ≤ bigM * (1 - m[:cn1_π_e_cc][c, q, eq[q].cc_id[i]]))
                @constraint(m, m[:p_e][c, q, bus] - eq[q].p[i] - (m[:σ][c] - 1) * eq[q].gen_slope[i] ≥ -bigM * (1 - m[:cn1_π_e_cc][c, q, eq[q].cc_id[i]]))
                @constraint(m, m[:p_e][c, q, bus] ≤ bigM * m[:cn1_π_e_cc][c, q, eq[q].cc_id[i]])
                @constraint(m, m[:p_e][c, q, bus] ≥ -bigM * m[:cn1_π_e_cc][c, q, eq[q].cc_id[i]])
                for (i, cc) in enumerate(eq[q].cc)
                    @constraint(m, m[:l_cc][c, q, i] == cc.l * m[:cn1_π_e_cc][c, q, i])
                end
            end
        end
    end
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

function loadloss!(m::Model, r::TNR, fixed_buses::Vector{VLabel}, eq::Vector{Equivalent})
    if isempty(eq)
        @constraint(m, [c in n_1cases(r)],
            m[:lostload][c] == sum(r.g[bus] - m[:load][c, bus] for bus in buses(r) if !(bus in fixed_buses) && r.g[bus] ≥ 0))
    else
        @constraint(m, [c in n_1cases(r)],
            m[:lostload][c] ==
            sum(r.g[bus] - m[:load][c, bus] for bus in buses(r) if !(bus in fixed_buses) && r.g[bus] ≥ 0) +
            sum(eq[q].cc[cc_id].l - m[:l_cc][c, q, cc_id] for q in eachindex(eq), cc_id in eachindex(eq[q].cc)))
        # sum(eq.p[i] - m[:p_e][c, bus] for (i, bus) in enumerate(eq.buses) if eq.p[i] ≥ 0))
    end
end

function branch_status!(m, r::TNR, branch_status::Union{AbstractArray{ELabel},Nothing})
    isnothing(branch_status) && return
    @constraint(m, [e in branches(r)], m[:v_branch][e...] == (e in branch_status))
end

function secured_dc_OTS(g::MetaGraph;
    contingencies::AbstractArray{ELabel}=ELabel[],
    n_1_connectedness::Bool=false,
    bus_confs::Vector{BusConf}=BusConf[],
    allow_branch_openings::Bool=true,
    OTS_only::Bool=false,
    tnr_pf::TNR_PF_TYPE=tnr_pf_phase,
    opf::Bool=false,
    bus_orig::String,
    fixed_buses::Vector{VLabel}=VLabel[],
    fixed_phases::Vector{Float64}=Float64[],
    branch_status::Union{AbstractArray{ELabel},Nothing}=nothing,
    equivalent::Union{Equivalent,Vector{Equivalent},Nothing}=nothing)

    eq = isnothing(equivalent) ? Equivalent[] :
         equivalent isa Equivalent ?
         Equivalent[equivalent] : equivalent
    model, r = init_model(g, contingencies, n_1_connectedness, bus_confs, allow_branch_openings, OTS_only, tnr_pf, opf, bus_orig)

    create_variables!(model, r, fixed_buses, eq)

    c_w!(model, r)

    if !r.opf
        N_balance!(model, r, fixed_buses, eq)
    end

    OTS_flows_phases!(model, r, fixed_buses, fixed_phases, eq)
    OTS_N_connectedness!(model, r, fixed_buses, eq)
    OTS_N_1_connectedness!(model, r, fixed_buses, eq)

    if r.opf
        OTS_balance_OPF!(model, r)
    else
        OTS_balance!(model, r, fixed_buses, eq)
    end

    loadloss!(model, r, fixed_buses, eq)

    overload!(model, r, eq)

    branch_status!(model, r, branch_status)

    @objective(model, Min,
        sum(model[:lostload][c] for c in n_1cases(r)) +
        00 * (r.opf ? sum(model[:loadshed][bus] for bus in buses(r)) : 0))
    optimize!(model)
    !is_solved_and_feasible(model) && @warn "NOT FEASIBLE"
    model, r
end
