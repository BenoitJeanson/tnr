include("equivalent.jl")
include("playUtils.jl")
include("equivalentTnrOptim.jl")
include("../src/tnr.jl")

# g, bus_confs, coord = create_mini_case();
# h = subgraph(g, Dict(2 => [1, 4], 4 => [2, 4]))

labs = collect('A':'Z')
contingencies = 1:20
trips = [2, 5, 4, 16, 10, 3, 6, 11, 19, 20, 1, 14, 15, 13, 12, 17, 18, 7, 9, 8] # used for griddraw
if !@isdefined g
    g, bus_confs, coord = create_case("case14", num -> "$(labs[num])")
    _, ϕ_base, _ = dc_flow_optim_fixed_phases!(g, [code_for(g, "A")], [0.0])
end

dsplit = Dict("E" => [("E", "F")], "D" => [("D", "I"), ("D", "G")])
fixed_bus = ["D", "E"]
bus_orig = "E"
sub, remaining = subgraph(g, dsplit)

# g, bus_confs, coord = create_mini_case()
# _, ϕ_base, _ = dc_flow_optim_fixed_phases!(g, [code_for(g, "1")], [0.0])
# bus_orig = "1"
# dsplit = Dict("2"=>[("1","2")], "4"=>[("1","4")])
# fixed_bus =["2", "4"]

# h=g
# bus_orig="A"

remaining_trip = [("A", "B"), ("A", "E"), ("B", "C"), ("B", "D"), ("B", "E"), ("C", "D"), ("D", "E")]
remaining_trip_ids = e_code_for(remaining, remaining_trip)

model, r = secured_dc_OTS(g,
    contingencies=1:ne(g),
    OTS_only=true,
    opf=false,
    bus_orig="A",
    # branch_status = convert_edgeids(result_g, g, h)
);
draw_TNR_result(g, model, r, coord=coord)



model, r = secured_dc_OTS(sub,
    contingencies=1:ne(sub),
    OTS_only=true,
    bus_confs=bus_confs,
    allow_branch_openings=true,
    tnr_pf=tnr_pf_phase,
    opf=false,
    bus_orig=bus_orig,
    fixed_buses=fixed_bus,
    fixed_phases=extract_export_phases(g, ϕ_base, fixed_bus),
    # branch_status = convert_edgeids(result_g, g, h)
)


# draw(h, layout=coord)



# flows, ϕ_base2, injections = dc_flow_optim_fixed_phases!(h, ["D", "E"], extract_export_phases(g, ϕ_base, ["D", "E"]))

equivalent = surrogate_parameters(sub, remaining, ["D", "E"], outages=[("J", "K"), ("I", "N"), ("G", "I")])
# eq = surrogate_parameters(sub, remaining, ["D", "E"])
# e_buses = convert_vertexids(to_code(g, fixed_bus), g, sub)


model, r = secured_dc_OTS(remaining,
    contingencies=1:ne(remaining),
    OTS_only=true,
    bus_confs=bus_confs,
    allow_branch_openings=true,
    tnr_pf=tnr_pf_phase,
    opf=false,
    bus_orig="A",
    # fixed_buses=fixed_bus,
    # fixed_phases=extract_export_phases(g, ϕ_base, fixed_bus),
    # branch_status = convert_edgeids(result_g, g, h)
    equivalent=equivalent
);