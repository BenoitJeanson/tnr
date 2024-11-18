using DrWatson
@quickactivate("tnr")

using Gurobi
env = Gurobi.Env()

include(joinpath(srcdir(), "DrawGrid.jl"))
include(joinpath(srcdir(), "GraphUtils.jl"))
include(joinpath(srcdir(), "dcpf.jl"))
include(joinpath(srcdir(), "tnr.jl"))
include(joinpath(srcdir(), "tnrOptim.jl"))

include("playUtils.jl")

g, bus_confs, coord = create_case("case14");

model, r = secured_dc_OTS(g,
                contingencies = 1:20,
                # contingencies = [1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16, 18, 19, 20],
                is_single_ρ = false,
                ρ_min_bound = 0.,
                n_1_connectedness = false,
                bus_confs = bus_confs,
                # bus_confs=BusConf[],
                allow_branch_openings = true,
                OTS_only = false
                )

optimize!(model);

trip = nothing #12
g_base = dc_flow(g, trip = trip)
layout = isa(coord, Dict{String, Point2}) ? coord : Stress(Ptype=Float32)
draw(g_base, trip = trip, layout = layout);

openings = Int[]
applied_conf=BusConf[]

report = []

if is_solved_and_feasible(model)
    store_result(model, "store_result.csv")
    push!(report, "overload: $(value(model[:overload]))")

    if r.allow_branch_openings || r.OTS_only
        openings = [i for (i, v_branch) in enumerate(value.(model[:v_branch])) if v_branch == 1]
        push!(report,"open: $openings")
        push!(report,[collect(edge_labels(g))[i] for i in openings])
    end
    
    if r.bus_confs ≠ BusConf[] && !r.OTS_only
        push!(report,"v_bus: $(value.(model[:v_bus]))")
        push!(report,"γ_bus: $(value.(model[:γ_bus]))")
        applied_conf = BusConf[r.bus_confs[i] for (i, v) in enumerate(value.(model[:v_bus])) if isapprox(v, 1, atol=1e-9)]
        push!(report,"applied_conf: $applied_conf")
    end
    # openings = [10]
    g_result, openings_result = add_subBus(g, applied_conf, openings)
    push!(report,"OOOOOOOOOPPPPPPPPPEEEEEEEEENNNNNNNIIIIIIIIINNNNNNNNNGGGGGGGG: $openings_result")
    dc_flow!(g_result, outages = [openings_result])
    draw(g_result , outages = openings_result, layout = layout);
else
    push!(report,"NOT FEASIBLE")
end

# show(stdout, MIME"text/plain",join(report,"\n"))
display(solution_summary(model))

# is_solved_and_feasible(model) && display(["$k: $(value(k))" for k in all_variables(model)])



# g, bus_confs, coord = create_case("case14");
# dc_flow!(g)
# fig = draw(g, layout = coord)
# save(joinpath("_research","tmp", "img.png"), fig)
# save(joinpath("_research","tmp", "img.svg"), fig)