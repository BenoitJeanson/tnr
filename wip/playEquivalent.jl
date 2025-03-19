include("../src/dcpf.jl")
include("Equivalent.jl")
include("EquivalentUtils.jl")
include("playUtils.jl")
include("equivalentTnrOptim.jl")
include("../src/tnr.jl")

labs = collect('A':'Z')
contingencies = 1:20
trips = [2, 5, 4, 16, 10, 3, 6, 11, 19, 20, 1, 14, 15, 13, 12, 17, 18, 7, 9, 8] # used for griddraw
bus_orig="A"
if !@isdefined g
    g, bus_confs, coord = create_case("case14", num -> "$(labs[num])")
end
_, ϕ_base, _ = dc_flow_optim!(g; fixed_buses=[code_for(g, "A")], fixed_phases=[0.0])

model_full, r_full = secured_dc_OTS(g,
    contingencies=1:ne(g),
    OTS_only=true,
    opf=false,
    bus_orig=bus_orig,
    # branch_status = convert_edgeids(result_g, g, h)
);
draw_TNR_result(g, model_full, r_full, coord=coord)

dsplit = Dict("E" => [("E", "F")], "D" => [("D", "I"), ("D", "G")])
fixed_buses = ["D", "E"]
bus_orig = "A"
sub_bus_orig = "E"
remaining_bus_orig = "A"
sub, remaining = subgraph(g, dsplit)

model_sub, r_sub = secured_dc_OTS(sub,
    contingencies=1:ne(sub),
    OTS_only=true,
    bus_orig=sub_bus_orig,
    fixed_buses=fixed_buses,
    fixed_phases=extract_export_phases(g, ϕ_base, fixed_buses),
);

eq_outage = equivalent_parameters(sub, remaining, fixed_buses, outages=[("F", "K"), ("I", "N"), ("G", "I")])
eq_loop = equivalent_parameters(sub, remaining, fixed_buses)

r = TNR(
    remaining,
    1:ne(remaining),
    # Int[1],
    false,
    BusConf[],
    true,
    true,
    tnr_pf_phase,
    false,
    remaining_bus_orig)

function sandbox(r, equivalent)
    m, r = init_model(r)
    _fixed_buses = Int[r.bus_orig_id]
    _fixed_phases = Float64[0.0]
    create_variables!(m, r, _fixed_buses, equivalent)
    c_w!(m, r)

    OTS_flows_phases!(m, r, _fixed_buses, _fixed_phases, equivalent)
    OTS_N_connectedness!(m, r, _fixed_buses, equivalent)

    OTS_N_1_connectedness!(m, r, _fixed_buses, equivalent)

    N_balance!(m, r, _fixed_buses, equivalent)
    OTS_balance!(m, r, _fixed_buses, equivalent)

    loadloss!(m, r, _fixed_buses, equivalent)

    overload!(m, r)

    # branch_status!(m, branch_status)

    @objective(m, Min,
        sum(m[:lostload][c] for c in n_1cases(r)))


    # @constraint(m, m[:v_branch_e][2] == 1)
    # @constraint(m, m[:v_branch][3] == 1)


    optimize!(m)
    if is_solved_and_feasible(m)
        @info "Model OK"
    else
        @warn "Model Infeasible"
    end
    m, r
end

m, _ = sandbox(r, eq_outage);
# m, _ = sandbox(r, eq_loop);
# m, _ = sandbox(r, nothing);


draw_TNR_result(remaining, m, r, coord=coord)
