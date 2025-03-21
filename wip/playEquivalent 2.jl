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

my_OTS(g, bus_orig; contingencies=1:ne(g), kwargs...) = secured_dc_OTS(g, contingencies=contingencies, OTS_only=true, bus_orig=bus_orig; kwargs...)

model_full, r_full = my_OTS(g, bus_orig);
draw_TNR_result(g, model_full, r_full; coord=coord)

dsplit2 = Dict("D" => [("D", "I")], "I" => [("G", "I"), ("D", "I")])
sub2, remaining2 = subgraph(g, dsplit2)
eq_buses = collect(keys(dsplit2))
sub_bus_orig2, remaining_bus_orig2 = "D", "A"

eq_sub = equivalent_parameters(sub2, remaining2, eq_buses, outages=[("G", "I")], true)
rem = copy(remaining2)
outages_rem = [("B", "E"), ("F", "K"), ("I", "N"), ("C", "D")]
calcanddraw(rem, outages=outages_rem, eq=eq_sub, layout=coord)
dc_flow_optim!(rem, outages=e_code_for(rem, outages_rem), equivalent=eq_sub,)
draw(rem, layout=coord, outages=e_code_for(rem, outages_rem))

eq_rem = equivalent_parameters(remaining2, sub2, eq_buses, outages=outages_rem, true)
sub = copy(sub2)
outages_sub = [("G", "I"), ("G", "H")]
calcanddraw(sub, outages=outages_sub, eq=eq_rem, layout=coord)
dc_flow_optim!(sub, outages=e_code_for(sub, outages_sub), equivalent=eq_rem,)
draw(sub, layout=coord, outages=e_code_for(sub, outages_sub))

outages_full = [outages_rem; outages_sub]
g_outages = copy(g)
dc_flow_optim!(g_outages, outages=e_code_for(g, outages_full))
draw(g_outages, layout=coord, outages=e_code_for(g_outages, outages_full))



# # _, ϕ_base, _ = dc_flow_optim!(g; fixed_buses=[code_for(g, "A")], fixed_phases=[0.0])
# # model_sub_fixed_ϕ, r_sub2_fixed_ϕ = my_OTS(sub2, sub_bus_orig2,
# #     fixed_buses=eq_buses, fixed_phases=extract_export_phases(g, ϕ_base, eq_buses));
# # draw_TNR_result(sub2, model_sub_fixed_ϕ, r_sub2_fixed_ϕ, coord=coord)

# ################################
# # eq_remaining2 = equivalent_parameters(remaining2, sub2, eq_buses, false)

# # model_sub2, r_sub2 = my_OTS(sub2, sub_bus_orig2, eq=eq_remaining2);
# # draw_TNR_result(sub2, model_sub2, r_sub2, coord=coord)
# # res_sub2 = e_label_for(sub2, openings(model_sub2))

# # eq_outage2 = equivalent_parameters(sub2, remaining2, eq_buses, outages=res_sub2, true)

# # model_rem2, r_rem2 = my_OTS(remaining2, remaining_bus_orig2, eq=eq_outage2);
# # draw_TNR_result(remaining2, model_rem2, r_rem2, coord=coord)
# # res_rem2 = e_label_for(remaining2, openings(model_rem2))

# # final_outages = [res_sub2; res_rem2]
# # griddraw(g, "A", trips, final_outages; layout=coord)

# # m_check, r_check = my_OTS(g, bus_orig, branch_status=final_outages);
# # draw_TNR_result(g, m_check, r_check, coord=coord)

# ################################
# eq_sub3 = equivalent_parameters(sub2, remaining2, eq_buses, true)

# model_remaining3, r_sub3 = my_OTS(remaining2, remaining_bus_orig2, eq=eq_sub3);
# draw_TNR_result(remaining2, model_remaining3, r_sub3; coord=coord)
# res_remaining3 = e_label_for(remaining2, openings(model_remaining3))

# eq_outage3 = equivalent_parameters(remaining2, sub2, eq_buses, outages=res_remaining3, true)

# model_sub3, r_sub3 = my_OTS(sub2, sub_bus_orig2, eq=eq_outage3);
# draw_TNR_result(sub2, model_sub3, r_sub3, coord=coord)
# res_sub3 = e_label_for(sub2, openings(model_sub3))

# final_outages = [res_remaining3; res_sub3]
# griddraw(g, "A", trips, final_outages; layout=coord)

# m_check, r_check = my_OTS(g, bus_orig, branch_status=final_outages);
# draw_TNR_result(g, m_check, r_check, coord=coord)



# model_sub_sand, r_sub_sand = sandbox(r_sub3, eq_outage3)
