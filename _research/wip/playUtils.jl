include(srcdir("DrawGrid.jl"))
include(srcdir("GraphUtils.jl"))
include(srcdir("dcpf.jl"))

function create_case(case::String, buslabels=lab -> "$lab")
    sys = System(joinpath("data", "exp_raw", "case", "$case.m"))
    coord = load_coord(joinpath("data", "exp_raw", "coord", "$case.csv"), buslabels)
    g = network2graph(sys, buslabels)
    bus_confs = Int[]

    if case == "case14"
        add_constraint(g, b -> b.p_max = 1.1)
        g[buslabels(14)] = 1
        g[buslabels(8)] = 0.1
        g["H"] = 0.5
        g["A", "B"].p_max = 3
        g["B", "D"].p_max = 3
        g["A", "E"].p_max = 3
        g["D", "E"].p_max = 1.7
        g["B", "E"].p_max = 0.9
        g["D", "I"].p_max = 0.33
        bus_confs = [BusConf(6, [SubBus(-0.17, [7, 9])])]

    elseif case == "case30"
        add_constraint(g, b -> b.p_max = .25)
        g["12","13"].p_max = .4
        g["6","8"].p_max = .4
        g["4","6"].p_max = .4
    end


    balance!(g, all_non_zero_uniform)

    g, bus_confs, coord
end

function create_mini_case(micro=false)
    g = build_simple_grid(micro=micro)
    g["1", "4"].p_max = 3
    # g["1","4"].p_max = 5
    # g["2"] = 2
    # g["1"] = -4

    mini_conf = [
        BusConf(3, [SubBus(0.5, [3])]),
        BusConf(2, [SubBus(0.5, [1, 3])]),
        BusConf(2, [SubBus(2, [1])])]

    g, mini_conf, load_coord(joinpath("data", "exp_raw", "coord", "mini.csv"))
end

function store_result(model, filename)
    open(filename, "w") do file
        foreach(k -> println(file, "$k, $(value(k))"), all_variables(model))
    end
end

function identifycriticalbranch(model)
    if !is_solved_and_feasible(model)
        @info "model not solved"
        return
    end
    ol = value(model[:overload]) + model.ext[:r].ρ_min_bound
    for i in eachindex(model[:flows])
        flow = model[:flows][i]
        if isapprox(value(flow), ol, atol=1e-8)
            return i
        end
    end
    return nothing
end

function displayflows(model)
    if !is_solved_and_feasible(model)
        @info "model not solved"
    end
    for flow in model[:flows]
        println("$flow: $(value(flow))")
    end
end

function showallvalues(model)
    if !is_solved_and_feasible(model)
        @info "model not solved"
        return
    end
    foreach(k -> println("$k: $(value(k))"), all_variables(model))
end

function calcanddraw(g; bus_orig=nothing, outages=Int[], trip=nothing, label="", kwargs...)
    bus_orig = isnothing(bus_orig) ? label_for(g, 1) : bus_orig
    thetrip = isa(trip, Tuple) ? e_code_for(g, trip...) : trip
    g_cc = connected(g, bus_orig, outages=outages, trip=thetrip)

    if nv(g) == nv(g_cc)
        h = g
        houtages = outages
        htrip = thetrip
    else
        h = g_cc
        balance!(h)
        houtages = Int[]
        htrip = nothing
        label = label == "" ? "" : L"\textbf{%$label}(%$(round(Int, 100*(total_load(g)-total_load(g_cc)))))"
    end

    dc_draw = dc_flow(h, outages=houtages, trip=htrip, pf_type=pf_linalg)
    draw(dc_draw, outages=houtages, trip=htrip, title=label,
        edge_width=br -> br.outage || br.trip ? 6 : 2,
        ; kwargs...)
end

function griddraw(g, bus_orig::String, trips::AbstractArray{Int}, outages=Int[]; kwargs...)
    fig = Figure(size=(600, 500), fontsize=20)
    splitlength = 5
    all_edges = collect(edges(g))
    for (i, trip) in enumerate(trips)
        e = all_edges[trip]
        label = "$(label_for(g, src(e)))$(label_for(g, dst(e)))"

        if !(trip in outages)
            draw(
                g;
                outages=outages,
                trip=trip,
                fig=fig[1, 1][(i-1)÷splitlength, (i-1)%splitlength],
                node_labels=nothing,
                edge_labels=nothing,
                arrow_size=0,
                node_size=nothing,
                edge_width=br -> br.trip || (abs(br.p > br.p_max + 1e-6)) ? 4 : 2,
                edge_coloring=br -> br.trip ? :black : br.outage ? :white : :black,
                kwargs...
                )

            calcanddraw(
                g;
                bus_orig=bus_orig,
                outages=outages,
                trip=trip,
                fig=fig[1, 1][(i-1)÷splitlength, (i-1)%splitlength],
                node_labels=nothing,
                edge_labels=nothing,
                arrow_size=0,
                node_size=nothing,
                edge_width=br -> br.trip || (abs(br.p > br.p_max + 1e-6)) ? 4 : 2,
                edge_coloring=br -> br.trip ? :black : br.outage ? :transparent : (abs(br.p > br.p_max + 1e-6)) ? :red : :green3,
                label=L"%$label", kwargs...)

        end
    end
    fig
end

function c_model_to_graph(g, model, orig_bus)
    if !is_solved_and_feasible(model)
        @info "model not fesible"
        return
    else
        h=deepcopy(g)
        openings = [i for (i, v_branch) in enumerate(value.(model[:v_branch])) if v_branch == 1]
        _edges = collect(edges(h))
        for (i, f) in enumerate(value.(model[:c_flows]))
            br = e_index_for(h, _edges[i])
            br.p = f/100
        end
        for v in vertices(h)
            h[label_for(g, v)] = 1/100
        end
        h[orig_bus] = -(nv(h)-1)/100

        draw(h , outages = openings, layout = layout, title = "Connectivity flows");
    end
end