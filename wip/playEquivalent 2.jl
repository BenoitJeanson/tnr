include("Equivalent.jl")
include("EquivalentUtils.jl")
include("playUtils.jl")
include("equivalentTnrOptim.jl")
include("../src/tnr.jl")

labs = collect('A':'Z')
contingencies = 1:20
bus_orig = "A"
trips = [2, 5, 4, 16, 10, 3, 6, 11, 19, 20, 1, 14, 15, 13, 12, 17, 18, 7, 9, 8] # used for griddraw
if !@isdefined g
    g, bus_confs, coord = create_case("case14", num -> "$(labs[num])")
end
_, ϕ_base, _ = dc_flow_optim!(g; fixed_buses=[code_for(g, "A")], fixed_phases=[0.0])

model_full, r_full = secured_dc_OTS(g,
    contingencies=1:ne(g),
    OTS_only=true,
    bus_orig=bus_orig,
);
draw_TNR_result(g, model_full, r_full, coord=coord)

dsplit2 = Dict("D" => [("D", "I")], "I" => [("G", "I"), ("D", "I")])
fixed_buses2 = ["D", "I"]
sub_bus_orig2 = "D"
remaining_bus_orig2 = "A"
sub2, remaining2 = subgraph(g, dsplit2)

model_sub_fixed_ϕ, r_sub2_fixed_ϕ = secured_dc_OTS(sub2,
    contingencies=1:ne(sub2),
    OTS_only=true,
    bus_orig=sub_bus_orig2,
    fixed_buses=fixed_buses2,
    fixed_phases=extract_export_phases(g, ϕ_base, fixed_buses2),
);

draw_TNR_result(sub2, model_sub_fixed_ϕ, r_sub2_fixed_ϕ, coord=coord)

eq_remaining2 = equivalent_parameters(remaining2, sub2, fixed_buses2, false)

model_sub2, r_sub2 = secured_dc_OTS(sub2,
    contingencies=1:ne(sub2),
    OTS_only=true,
    bus_orig=sub_bus_orig2,
    eq=eq_remaining2
);

draw_TNR_result(sub2, model_sub2, r_sub2, coord=coord)

eq_outage2 = equivalent_parameters(sub2, remaining2, fixed_buses2, outages=[("D", "I")], true)
# eq_loop = equivalent_parameters(sub, remaining, fixed_buses)

model_rem2, r_rem2 = secured_dc_OTS(remaining2,
    contingencies=1:ne(remaining2),
    OTS_only=true,
    bus_orig=remaining_bus_orig2,
    eq=eq_outage2
);

draw_TNR_result(remaining2, model_rem2, r_rem2, coord=coord)

final_outages = [("B", "E"), ("I", "N"), ("I", "J"), ("C", "D"), ("D", "I")]
griddraw(g, "A", trips, final_outages; layout=coord)

m_check, r_check = secured_dc_OTS(g, contingencies=1:ne(g), OTS_only=true, bus_orig=bus_orig,
    branch_status=[e_label_for(g, e) in final_outages for e in 1:ne(g)]
)

draw_TNR_result(g, m_check, r_check, coord=coord)
