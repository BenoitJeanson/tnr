include("Equivalent.jl")
include("EquivalentUtils.jl")
include("playUtils.jl")
include("equivalentTnrOptim.jl")
include("../src/tnr.jl")

labs = collect('A':'Z')
bus_orig = "A"
if !@isdefined g
    g, trippings, bus_confs, coord = create_case("case14", num -> "$(labs[num])")
end

my_OTS(g, bus_orig; contingencies=collect(edge_labels(g)), kwargs...) = secured_dc_OTS(g, contingencies=contingencies, OTS_only=true, bus_orig=bus_orig; kwargs...)

model_full, r_full = my_OTS(g, bus_orig);
draw_TNR_result(g, model_full, r_full; coord=coord)

dsplit = Dict("D" => [("D", "I")], "I" => [("G", "I"), ("D", "I")])
sub, remaining = subgraph(g, dsplit)
buses_eq = collect(keys(dsplit))
sub_bus_orig, remaining_bus_orig = "D", "A"

eq_sub3 = equivalent_parameters(sub, buses_eq, true)

model_remaining3, r_sub3 = my_OTS(remaining, remaining_bus_orig, eq=eq_sub3);
draw_TNR_result(remaining, model_remaining3, r_sub3; coord=coord)
res_remaining3 = openings(model_remaining3)

eq_outage3 = equivalent_parameters(remaining, buses_eq, outages=res_remaining3, true)

model_sub3, r_sub3 = my_OTS(sub, sub_bus_orig, eq=eq_outage3);
draw_TNR_result(sub, model_sub3, r_sub3, coord=coord)
res_sub3 = openings(model_sub3)

final_outages = [res_remaining3; res_sub3]
griddraw(g, "A", trippings, final_outages; layout=coord)

m_check, r_check = my_OTS(g, bus_orig, branch_status=final_outages);
draw_TNR_result(g, m_check, r_check, coord=coord)

