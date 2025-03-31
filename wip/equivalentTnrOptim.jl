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

function create_variables!(m::Model, r::TNR, fixed_buses::Vector{VLabel}=VLabel[], eq::Union{Equivalent,Nothing}=nothing)

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

    if !isnothing(eq)
        @variable(m, flows_e[cases(r), eq.buses])
        @variable(m, inner_flows_e[cases(r), [str(e) for e in eq.internal_branches]])
        if eq.allow_connector_opening
            @variable(m, ϕ_e[cases(r), eq.buses])
            @variable(m, v_branch_e[eq.buses], Bin)
            @variable(m, cn1_ψ_e_bus[n_1cases(r), e_bus in eq.buses], Bin)
            @variable(m, cn1_ψ_e_cc[n_1cases(r), e_cc in 1:length(eq.cc), eq.cc[e_cc].buses], Bin)
        end
        @variable(m, p_e[cases(r), eq.buses])
        @variable(m, l_cc[n_1cases(r), 1:length(eq.cc)])


        @variable(m, c_flows_e[r.n_1_connectedness ? cases(r) : [n_case(r)], eq.buses])

        @variable(m, cn1_flows_e[n_1cases(r), eq.buses])
        @variable(m, cn1_π_e_cc[n_1cases(r), 1:length(eq.cc)], Bin)
    end

    @variable(m, v_branch[busfrom in labels(r.g), outneighbor_labels(r.g, busfrom)], Bin)

    @variable(m, ϕ[cases(r), buses(r)])

    if r.opf
        @variable(m, loadshed[buses(r)] ≥ 0)
        @variable(m, σ[cases(r)] ≥ 0)
    else
        @variable(m, σ[n_1cases(r)] ≥ 0)
    end

end

function OTS_flows_phases!(m::Model, r::TNR, fixed_buses::Vector{VLabel}, fixed_ϕ::Vector{Float64}, eq::Union{Equivalent,Nothing})
    big_M = 100 #TODO: to change  
    @info "HARDCODED OTS_flows_phases bigM to $big_M"

    @constraint(m, [c in cases(r), (i, f_bus) in enumerate(fixed_buses)], m[:ϕ][c, f_bus] == fixed_ϕ[i])

    @constraint(m, [c in cases(r)], m[:flows][c, :, :] .≤ (1 .- m[:c_w][c, :, :]) .* big_M)
    @constraint(m, [c in cases(r)], -m[:flows][c, :, :] .≤ (1 .- m[:c_w][c, :, :]) .* big_M)

    @constraint(m, [c in cases(r), b in branches(r)],
        m[:flows][c, b...] - sum(r.g[b...].b * (m[:ϕ][c, to(b)] - m[:ϕ][c, from(b)])) ≤ big_M * m[:c_w][c, b...])
    @constraint(m, [c in cases(r), b in branches(r)],
        m[:flows][c, b...] - sum(r.g[b...].b * (m[:ϕ][c, to(b)] - m[:ϕ][c, from(b)])) ≥ -big_M * m[:c_w][c, b...])

    if isnothing(eq)
        @constraint(m, [c in cases(r), bus in buses(r); !(bus in fixed_buses)],
            m[:load][c, bus] - m[:gen][c, bus] ==
            sum(sign * m[:flows][c, b...] for (b, sign) in incident_signed(r, bus)))
    else
        @constraint(m, [c in cases(r), bus in buses(r); !(bus in fixed_buses)],
            m[:load][c, bus] - m[:gen][c, bus] ==
            sum(sign * m[:flows][c, b...] for (b, sign) in incident_signed(r, bus)) +
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
                    m[:flows_e][c, bus] - m[:p_e][c, bus] ==
                    sum(eq.B[bus_id, bus_ϕ_id] * m[:ϕ][c, bus_ϕ] for (bus_ϕ_id, bus_ϕ) in enumerate(eq.buses)))
            end
        end
        @constraint(m, [c in cases(r), cc in eq.cc],
            sum(m[:flows_e][c, bus] for bus in cc.buses) ==
            sum(m[:p_e][c, bus] for bus in cc.buses))
    end
end

function overload!(m::Model, r::TNR, eq::Union{Equivalent,Nothing})
    @constraint(m, [c in cases(r), b in branches(r)], m[:flows][c, b...] .≤ r.g[b...].p_max)
    @constraint(m, [c in cases(r), b in branches(r)], -m[:flows][c, b...] .≤ r.g[b...].p_max)

    if !isnothing(eq)
        if eq.allow_connector_opening
            @constraint(m, [c in cases(r), e_id in 1:length(eq.f)], eq.f[e_id] + sum(eq.F[e_id, :] .* m[:ϕ_e][c, :]) ≤ eq.f_max[e_id])
            @constraint(m, [c in cases(r), e_id in 1:length(eq.f)], eq.f[e_id] + sum(eq.F[e_id, :] .* m[:ϕ_e][c, :]) ≥ -eq.f_max[e_id])
        else
            ib = [str(eq.internal_branches[e_id]) for e_id in 1:length(eq.f)]
            @constraint(m, [c in cases(r), e_id in 1:length(eq.f)],
                m[:inner_flows_e][c, ib[e_id]] == eq.f[e_id] + sum(eq.F[e_id, bus_id] * m[:ϕ][c, bus] for (bus_id, bus) in enumerate(eq.buses)))
            @constraint(m, [c in cases(r), e_id in 1:length(eq.f)],
                m[:inner_flows_e][c, ib[e_id]] ≤ eq.f_max[e_id])
            @constraint(m, [c in cases(r), e_id in 1:length(eq.f), bus in eq.buses],
                m[:inner_flows_e][c, ib[e_id]] ≥ -eq.f_max[e_id])
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

function OTS_N_connectedness!(m::Model, r::TNR, fixed_buses::Vector{VLabel}, eq::Union{Nothing,Equivalent})
    _fixed_buses = isempty(fixed_buses) ? VLabel[r.bus_orig] : fixed_buses
    big_M = nb_buses(r) + (isnothing(eq) ? 0 : length(eq.buses))
    cases = r.n_1_connectedness ? cases(r) : [n_case(r)]

    @constraint(m, [c in cases], m[:c_flows][c, :, :] .≤ big_M .* (1 .- m[:c_w][c, :, :]))
    @constraint(m, [c in cases], -m[:c_flows][c, :, :] .≤ big_M .* (1 .- m[:c_w][c, :, :]))

    if isnothing(eq)
        @constraint(m, [c in cases, bus in buses(r); !(bus in _fixed_buses)],
            sum(sign * m[:c_flows][c, b...] for (b, sign) in incident_signed(r, bus)) == 1)
        @constraint(m, [c in cases, bus in _fixed_buses],
            sum(sign * m[:c_flows][c, b...] for (b, sign) in incident_signed(r, bus)) == m[:c_slacks][bus])
        @constraint(m, sum(m[:c_slacks][bus] for bus in _fixed_buses) == length(_fixed_buses) - nb_buses(r))

    else
        @constraint(m, [c in cases, bus in buses(r); !(bus in _fixed_buses)],
            sum(sign * m[:c_flows][c, b...] for (b, sign) in incident_signed(r, bus)) -
            ((bus in eq.buses) ? m[:c_flows_e][c, bus] : 0) == 1)
        @constraint(m, [c in cases, bus in _fixed_buses],
            sum(sign * m[:c_flows][c, b...] for (b, sign) in incident_signed(r, bus)) -
            ((bus in eq.buses) ? m[:c_flows_e][c, bus] : 0) == m[:c_slacks][bus])
        @constraint(m, sum(m[:c_slacks][bus] for bus in _fixed_buses) == length(_fixed_buses) - nb_buses(r) - length(eq.cc))

        if eq.allow_connector_opening
            @constraint(m, [c in cases, e_bus in eq.buses], m[:c_flows_e][c, e_bus] ≤ big_M * (1 - m[:v_branch_e][e_bus]))
            @constraint(m, [c in cases, e_bus in eq.buses], -m[:c_flows_e][c, e_bus] ≤ big_M * (1 - m[:v_branch_e][e_bus]))
        end

        @constraint(m, [c in cases, cc in eq.cc], sum(m[:c_flows_e][c, br] for br in cc.buses) == 1)
    end

end

function OTS_N_1_connectedness!(m::Model, r::TNR, fixed_buses::Vector{VLabel}, eq::Union{Nothing,Equivalent})
    bigM_nb_v = nb_buses(r) + 2 + (isnothing(eq) ? 0 : length(eq.buses))

    _fixed_buses = isempty(fixed_buses) ? VLabel[r.bus_orig] : fixed_buses

    # Set of constraints to ensure cn1_π is 0 in unenergized ares.
    @constraint(m, [c in n_1cases(r), b in branches(r)], m[:cn1_flows][c, b...] .≤ bigM_nb_v .* (1 .- m[:c_w][c, b...]))
    @constraint(m, [c in n_1cases(r), e in branches(r)], -m[:cn1_flows][c, :, :] .≤ bigM_nb_v .* (1 .- m[:c_w][c, :, :]))

    if isnothing(eq)
        @constraint(m, [c in n_1cases(r), bus in [bus for bus in buses(r) if !(bus in _fixed_buses)]],
            sum(sign * m[:cn1_flows][c, b...] for (b, sign) in incident_signed(r, bus)) ==
            m[:cn1_π][c, bus])
    else
        @constraint(m, [c in n_1cases(r), bus in [bus for bus in buses(r) if !(bus in _fixed_buses)]],
            sum(sign * m[:cn1_flows][c, b...] for (b, sign) in incident_signed(r, bus)) -
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

    @constraint(m, [c in n_1cases(r), bus in buses(r), b in incident(r, bus); !(bus in _fixed_buses)],
        m[:cn1_π][c, bus] ≥ m[:cn1_ψ][c, bus, b...])

    @constraint(m, [c in n_1cases(r), bus in buses(r), b in incident(r, bus); !(bus in _fixed_buses)],
        m[:cn1_ψ][c, bus, b...] ≥ m[:cn1_π][c, opposite(b, bus)] - m[:c_w][c, b...])
    @constraint(m, [c in n_1cases(r), bus in buses(r), b in incident(r, bus); !(bus in _fixed_buses)],
        m[:cn1_ψ][c, bus, b...] ≤ m[:cn1_π][c, opposite(b, bus)])
    @constraint(m, [c in n_1cases(r), bus in buses(r), b in incident(r, bus); !(bus in _fixed_buses)],
        m[:cn1_ψ][c, bus, b...] ≤ 1 - m[:c_w][c, b...])

    if isnothing(eq)
        @constraint(m, [c in n_1cases(r), bus in buses(r); !(bus in _fixed_buses)],
            m[:cn1_π][c, bus] ≤ sum(m[:cn1_ψ][c, bus, b...] for b in incident(r, bus)))
    else
        if eq.allow_connector_opening
            @constraint(m, [c in n_1cases(r), bus in buses(r); !(bus in _fixed_buses)],
                m[:cn1_π][c, bus] ≤
                sum(m[:cn1_ψ][c, bus, b...] for b in incident(r, bus)) +
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
                sum(m[:cn1_ψ][c, bus, br...] for br in incident(r, bus)) +
                (bus in eq.buses ? m[:cn1_π_e_cc][c, get_cci(eq, bus)] : 0))
            @constraint(m, [c in n_1cases(r), bus in eq.buses],
                m[:cn1_π][c, bus] == m[:cn1_π_e_cc][c, get_cci(eq, bus)])
        end
    end
end

function N_balance!(m::Model, r::TNR, fixed_buses::Vector{VLabel}, eq::Union{Nothing,Equivalent})
    for bus in buses(r)
        (bus in fixed_buses) && continue
        @constraint(m, m[:load][n_case(r), bus] == max(0, r.g[bus]))
        @constraint(m, m[:gen][n_case(r), bus] == -min(0, r.g[bus]))
    end
    if !isnothing(eq)
        for (i, bus) in enumerate(eq.buses)
            @constraint(m, m[:p_e][n_case(r), bus] == eq.p[i])
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

function OTS_balance!(m::Model, r::TNR, fixed_buses::Vector{VLabel}, eq::Union{Nothing,Equivalent})
    bigM = 10 # TODO: to define
    @info "HARDCODED OTS_balance - bigM to $bigM"

    if isempty(fixed_buses)
        if isnothing(eq)
            @constraint(m, [c in n_1cases(r)], sum(m[:load][c, :] .- m[:gen][c, :]) == 0)
        else
            @constraint(m, [c in n_1cases(r)], sum(m[:load][c, :] .- m[:gen][c, :]) - sum(m[:p_e][c, :]) == 0)
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
    end
    if !isnothing(eq)
        for c in n_1cases(r)
            for (i, bus) in enumerate(eq.buses)
                @constraint(m, m[:p_e][c, bus] - eq.p[i] - (m[:σ][c] - 1) * eq.gen_slope[i] ≤ bigM * (1 - m[:cn1_π_e_cc][c, eq.cc_id[i]]))
                @constraint(m, m[:p_e][c, bus] - eq.p[i] - (m[:σ][c] - 1) * eq.gen_slope[i] ≥ -bigM * (1 - m[:cn1_π_e_cc][c, eq.cc_id[i]]))
                @constraint(m, m[:p_e][c, bus] ≤ bigM * m[:cn1_π_e_cc][c, eq.cc_id[i]])
                @constraint(m, m[:p_e][c, bus] ≥ -bigM * m[:cn1_π_e_cc][c, eq.cc_id[i]])
            end
            for (i, cc) in enumerate(eq.cc)
                @constraint(m, m[:l_cc][c, i] == cc.l * m[:cn1_π_e_cc][c, i])
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

function loadloss!(m::Model, r::TNR, fixed_buses::Vector{VLabel}, eq::Union{Nothing,Equivalent})
    if isnothing(eq)
        @constraint(m, [c in n_1cases(r)],
            m[:lostload][c] == sum(r.g[bus] - m[:load][c, bus] for bus in buses(r) if !(bus in fixed_buses) && r.g[bus] ≥ 0))
    else
        @constraint(m, [c in n_1cases(r)],
            m[:lostload][c] ==
            sum(r.g[bus] - m[:load][c, bus] for bus in buses(r) if !(bus in fixed_buses) && r.g[bus] ≥ 0) +
            sum(eq.cc[cc_id].l - m[:l_cc][c, cc_id] for cc_id in 1:length(eq.cc)))
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
    eq::Union{Equivalent,Nothing}=nothing)

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
