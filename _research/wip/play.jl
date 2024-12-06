using DrWatson
@quickactivate("tnr")

using Gurobi
env = Gurobi.Env()

include(srcdir("DrawGrid.jl"))
include(srcdir("GraphUtils.jl"))
include(srcdir("dcpf.jl"))
include(srcdir("tnrOptim.jl"))

include("playUtils.jl")

g, bus_confs, coord = create_case("case14");

# g, bus_confs, coord = create_mini_case()

model, r = secured_dc_OTS(g,
                contingencies = 1:20,
                # contingencies = [1, 2, 14, 15],
                is_single_ρ = true,
                ρ_min_bound = 1.,
                # bus_confs = bus_confs,
                allow_branch_openings = true,
                OTS_only = true
                )

optimize!(model);

fig = Figure(size = (1500, 800), fontsize = 15)

trip = nothing #12
g_base = dc_flow(g, trip = trip)
layout = isa(coord, Dict{String, Point2}) ? coord : Stress(Ptype=Float32)
draw(g_base, fig = fig[1,1], trip = trip, layout = layout);

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
    draw(g_result , fig = fig[1,2], outages = openings_result, layout = layout);
    println("OK")
else
    println("NOT FEASIBLE")
end

display(fig)

# is_solved_and_feasible(model) && display(["$k: $(value(k))" for k in all_variables(model)])



# g, bus_confs, coord = create_case("case14");
# dc_flow!(g)
# fig = draw(g, layout = coord)
# save(projectdir("_research","tmp", "img.png"), fig)
# save(joinpath("_research","tmp", "img.svg"), fig)