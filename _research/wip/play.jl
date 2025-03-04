using DrWatson
@quickactivate("tnr")

using Gurobi
env = Gurobi.Env()

include(srcdir("DrawGrid.jl"))
include(srcdir("GraphUtils.jl"))
include(srcdir("dcpf.jl"))
include(srcdir("tnrOptim.jl"))

include("playUtils.jl")

labs = collect('A':'Z')
g, bus_confs, coord = create_case("case14", num -> "$(labs[num])");
contingencies = 1:20
trips = [2, 5, 4, 16, 10, 3, 6, 11, 19, 20, 1, 14, 15, 13, 12, 17, 18, 7, 9, 8]
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

optimize!(model);

fig = Figure(size = (1500, 1600), fontsize = 20)

trip = nothing #12
g_base = dc_flow(g, trip = trip)
layout = isa(coord, Dict{String, Point2}) ? coord : Stress(Ptype=Float32)
draw(g_base, fig = fig[1,1], trip = trip, layout = layout, title = "$(r.tnr_pf)");

openings = Int[]
applied_conf=BusConf[]

report = []

display(solution_summary(model))

if is_solved_and_feasible(model)
    store_result(model, "store_result.csv")
    println("Overload: $(value(model[:overload]))")

    if r.allow_branch_openings || r.OTS_only
        push!
        openings = [i for (i, v_branch) in enumerate(value.(model[:v_branch])) if v_branch == 1]
        println("Open: $openings")
        println([collect(edge_labels(g))[i] for i in openings])
    end
    
    if r.bus_confs ≠ BusConf[] && !r.OTS_only
        println("v_bus: $(value.(model[:v_bus]))")
        applied_conf = BusConf[r.bus_confs[i] for (i, v) in enumerate(value.(model[:v_bus])) if isapprox(v, 1, atol=1e-9)]
        println("Applied_conf: $applied_conf")
    end
    g_result, openings_result = add_subBus(g, applied_conf, openings)
    dc_flow!(g_result, outages = [openings_result])
    draw(g_result , fig = fig[1,2], outages = openings_result, layout = layout, title = "objective value:$(round(objective_value(model); digits = 5 ))");
    println("OK")
else
    println("NOT FEASIBLE")
end

display(fig)

# is_solved_and_feasible(model) && display(["$k: $(value(k))" for k in all_variables(model)])
