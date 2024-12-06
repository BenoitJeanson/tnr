function create_case(case::String)
    sys = System(joinpath("data","exp_raw","case", "$case.m"));
    coord = load_coord(joinpath("data","exp_raw","coord", "$case.csv"))
    g = network2graph(sys)
    g["14"] = 1
    balance!(g, all_non_zero_uniform)

    add_constraint(g, b->b.p_max=1.3)
    if case =="case14"
        g["1","2"].p_max = 2.3
        g["1","5"].p_max = 2.3
        # g["5","6"].p_max = 2.3
        # g["4","5"].p_max = 2.3
        # g["5","6"].p_max = .9
        # g["4","7"].p_max = .6
        # g["4","9"].p_max = .6
    end

    bus_confs = [BusConf(6, [SubBus(-.17 , [7, 9])])]
    g, bus_confs, coord
end

function store_result(model, filename)
    open(filename, "w") do file
        foreach(k->println(file, "$k, $(value(k))"), all_variables(model))
    end
end

function create_mini_case(micro=false)
    g = build_simple_grid(micro=micro)
    g["1","4"].p_max = 3
    # g["1","4"].p_max = 5
    # g["2"] = 2
    # g["1"] = -4

    mini_conf = [
        BusConf(3, [SubBus(.5, [3])]),
        BusConf(2, [SubBus(.5, [1, 3])]),
        BusConf(2, [SubBus(2, [1])])]
    
    g, mini_conf, load_coord(joinpath("data","exp_raw","coord", "mini.csv"))
end

function identifycriticalbranch(model)
    if !is_solved_and_feasible(model)
        println("model not solved")
        return
    end
    ol = value(model[:overload]) + model.ext[:r].ρ_min_bound
    for i in eachindex(model[:flows])
        flow = model[:flows][i]
        if isapprox(value(flow), ol, atol = 1e-8)
            return i
        end
    end
    return nothing
end

function displayflows(model)
    if !is_solved_and_feasible(model)
        println("model not solved")
    end
    for flow in model[:flows]
        println("$flow: $(value(flow))")
    end
end

function showallvalues(model)
    if !is_solved_and_feasible(model)
        println("model not solved")
        return
    end
    foreach(k->println("$k: $(value(k))"), all_variables(model))
end

function calcanddraw(g; outages=Int[], trip=nothing, kwargs...)
    dc_draw = dc_flow(g, outages = outages, trip=trip)
    draw(dc_draw, outages = outages, trip=trip; kwargs...) 
end

function griddraw(g, bus_orig::String, trips::AbstractArray{Int}, outages=Int[]; kwargs...)
    fig = Figure(size = (1500, 1500), fontsize = 30)
    calcanddraw(
        g,
        outages = outages,
        trip = trip,
        fig = fig[1,1],
        edge_labels = nothing,
        node_labels = (g, label; digits = 2)->label,
        title = "Base, bus origign = $bus_orig",
        layout = layout, kwargs...)
    
    splitlength = trunc(Int, sqrt(length(trips)))
    all_edges = collect(edges(g))
    for (i,trip) in enumerate(trips)
        g_cc = connected(g, bus_orig, outages=outages, trip=trip)
        e = all_edges[i]
        label = "$(label_for(g, src(e)))-$(label_for(g, dst(e)))"
        if nv(g) == nv(g_cc)
            h = g
            houtages = outages
            htrip = trip
        else
            h = g_cc
            balance!(h)
            houtages = Int[]
            htrip = nothing
            label = "$label *"
        end
        calcanddraw(
            h,
            outages = houtages,
            trip = htrip,
            fig = fig[2:3,1][(i-1)÷splitlength, (i-1)%splitlength],
            node_labels = nothing,
            edge_labels = nothing,
            arrow_size = 0,
            node_size = nothing,
            title = "$i: $label",
            layout = layout, kwargs...)
    end
    fig
end
