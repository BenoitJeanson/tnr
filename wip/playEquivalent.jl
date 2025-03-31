include("../src/dcpf.jl")
include("../src/tnr.jl")
include("playUtils.jl")
include("Equivalent.jl")
include("EquivalentUtils.jl")
include("equivalentTnrOptim.jl")

labs = collect('A':'Z')
contingencies = 1:20
bus_orig = "A"
if !@isdefined g
    g, trippings, bus_confs, coord = create_case("case14", num -> "$(labs[num])")
end
_, ϕ_base, _ = dc_flow_optim!(g; slack_buses=["A"], slack_phases=[0.0])

model_full, r_full = secured_dc_OTS(g,
    contingencies=collect(edge_labels(g)),
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
    contingencies=collect(edge_labels(sub)),
    OTS_only=true,
    bus_orig=sub_bus_orig,
    fixed_buses=fixed_buses,
    fixed_phases=extract_export_phases(g, ϕ_base, fixed_buses),
);
draw_TNR_result(sub, model_sub, r_sub, slack_buses=fixed_buses, slack_phases=extract_export_phases(g, ϕ_base, fixed_buses), coord=coord)

eq_rem = equivalent_parameters(remaining, "remaining", fixed_buses, true)
model_sub_eq, r_sub_eq = secured_dc_OTS(sub,
    contingencies=collect(edge_labels(sub)),
    OTS_only=true,
    bus_orig=sub_bus_orig,
    equivalent = eq_rem
);
draw_TNR_result(sub, model_sub_eq, r_sub_eq, coord=coord)

eq_outage = equivalent_parameters(sub, "sub_outage", fixed_buses, false, outages=openings(model_sub_eq))
eq_loop = equivalent_parameters(sub, "sub looped", fixed_buses, true)


m, r = secured_dc_OTS(remaining,
    contingencies=collect(edge_labels(remaining)),
    OTS_only=true,
    bus_orig="A",
    # fixed_buses=["A"],
    # fixed_phases=[0.],
    equivalent=eq_outage)
draw_TNR_result(remaining, m, r, coord=coord)
