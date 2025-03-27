using Gurobi
env = Gurobi.Env()

include("../src/DrawGrid.jl")
include("../src/GraphUtils.jl")
include("../src/dcpf.jl")
include("../src/tnrOptim.jl")

include("playUtils.jl")

labs = collect('A':'Z')
if !@isdefined g
    g, bus_confs, coord = create_case("case14", num -> "$(labs[num])");
end
contingencies = 1:20
trips = [2, 5, 4, 16, 10, 3, 6, 11, 19, 20, 1, 14, 15, 13, 12, 17, 18, 7, 9, 8] # used for griddraw
bus_orig = "A"

# g, bus_confs, coord = create_case("case30", num -> "$(num)");

# g, bus_confs, coord = create_mini_case()
# contingencies = 1:5
# bus_orig = "1"

# nb_case=4
# g, coord = create_multiplecase(nb_case)
# contingencies = 1:(20*nb_case)

tnr_pf = tnr_pf_phase
@info "Secured DC OTS, with $tnr_pf"
model, r = secured_dc_OTS(g,
                contingencies = contingencies,
                OTS_only = true,
                bus_confs = bus_confs,
                allow_branch_openings = true,
                tnr_pf = tnr_pf,
                opf = true,
                bus_orig = bus_orig
                )

draw_TNR_result(g, model, r, coord=coord)