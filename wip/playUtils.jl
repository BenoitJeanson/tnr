include("../src/DrawGrid.jl")
include("../src/GraphUtils.jl")
include("../src/dcpf.jl")

function create_case(case::String, buslabels=lab -> "$lab")
    sys = System(joinpath("data", "exp_raw", "case", "$case.m"))
    coord = load_coord(joinpath("data", "exp_raw", "coord", "$case.csv"), buslabels)
    g = network2graph(sys, buslabels)
    bus_confs = Int[]
    trippings=nothing

    if case == "case14"
        trippings = [("F", "L"),("L", "M"),("F", "M"),("M", "N"),("I", "N"),("F", "K"),("J", "K"),("I", "J"),("G", "I"),("G", "H"),("E", "F"),("D", "I"),("D", "G"),("D", "E"),("C", "D"),("A", "E"),("A", "B"),("B", "E"),("B", "D"),("B", "C")] # used for griddra
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
        add_constraint(g, b -> b.p_max = 0.25)
        g["12", "13"].p_max = 0.4
        g["6", "8"].p_max = 0.4
        g["4", "6"].p_max = 0.4
    end


    balance!(g, all_non_zero_uniform)

    g, trippings, bus_confs, coord
end

function create_multiplecase(n::Int=2)
    @info "create create_multiplecase, size = $n"

    labs = collect('A':'Z')
    g_base, _, coord = create_case("case14", num -> "$(labs[num])")
    g = copy(g_base)
    for i in 2:n
        for v in vertices(g_base)
            lab = label_for(g_base, v)
            if lab ≠ "A"
                g["$(lab)$i"] = g_base[lab]
                pt = coord[lab]
                coord["$(lab)$i"] = Point(i ≤ 2 ? pt[1] : (-pt[1] - 2), i == 3 ? pt[2] : -pt[2] - 2)
            end
        end
        for e in edges(g_base)
            labs, labd = lsrc(g_base, e), ldst(g_base, e)
            g[labs == "A" ? "A" : "$(labs)$i", labd == "A" ? "A" : "$(labd)$i"] = deepcopy(g_base[labs, labd])
        end
    end
    coord["A"] = Point(n > 2 ? -1 : 0, -1)
    g["A"] = n * g["A"]
    g, coord
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

function connected(g::MetaGraph, bus_origin::String; outages::AbstractArray{ELabel}=ELabel[], trip::Union{Nothing,ELabel}=nothing, eq::Union{Nothing,Equivalent}=nothing)

    """
    _extractsubgraph is embedded as it only applies if the bus_labels are the labels of buses belonging to a single connected component 
    """
    function _extractsubgraph(g::MetaGraph, bus_labels::Vector{String})
        g_sub = _initgraph()
        for bus in bus_labels
            g_sub[bus] = g[bus]
        end
        for e in edges(g)   # TODO: prefere use of edge_labels
            slabel = label_for(g, src(e))
            if slabel in bus_labels
                dlabel = label_for(g, dst(e))
                g_sub[slabel, dlabel] = g[slabel, dlabel]
            end
        end
        g_sub
    end

    g_result = copy(g)

    openbranches = copy(outages)
    !isnothing(trip) && push!(openbranches, trip)

    for br in openbranches
        delete!(g_result, br...)
    end

    # create virtual edges that will be removed after the connected components calculation
    eq_edges = ELabel[]
    if !isnothing(eq)
        for cc in eq.cc
            for i in eachindex(cc.buses), j in i+1:length(cc.buses)
                elabs = (cc.buses[i], cc.buses[j])
                if !haskey(g_result, elabs...)
                    g_result[elabs...] = Branch()
                    push!(eq_edges, elabs)
                end
            end
        end
    end

    ccs = connected_components(g_result)
    foreach(elab -> delete!(g_result, elab...), eq_edges)
    cc = ccs[findfirst(cc -> bus_origin in [label_for(g_result, v) for v in cc], ccs)]
    _extractsubgraph(g_result, [label_for(g_result, v) for v in cc])
end

function total_load(g)
    result = 0
    for v in vertices(g)
        lab = label_for(g, v)
        if g[lab] > 0
            result += g[lab]
        end
    end
    result
end


function calcanddraw(g; bus_orig=nothing,
    outages::AbstractArray{ELabel}=ELabel[],
    trip::Union{Nothing,ELabel}=nothing,
    label="", eq=nothing, kwargs...)

    bus_orig = isnothing(bus_orig) ? label_for(g, 1) : bus_orig
    g_cc = connected(g, bus_orig; outages=outages, trip=trip, eq=eq)

    if nv(g) == nv(g_cc)
        h = g
        houtages = outages
        htrip = trip
    else
        h = g_cc
        balance!(h)
        houtages = ELabel[]
        htrip = nothing
        label = label == "" ? "" : L"\textbf{%$label}(%$(round(Int, 100*(total_load(g)-total_load(g_cc)))))"
    end
    dc_draw = copy(h)
    dc_flow_optim!(dc_draw; outages=houtages, trip=htrip, equivalent=eq)
    draw(dc_draw, outages=houtages, trip=htrip, title=label,
        edge_width=br -> br.outage || br.trip ? 6 : 2,
        ; kwargs...)
end

function griddraw(g, bus_orig::String, trips::Union{AbstractArray{ELabel}}, outages::Union{AbstractArray{ELabel}}=ELabel[], eq=nothing; kwargs...)
    fig = Figure(size=(600, 500), fontsize=20)
    splitlength = 5
    for (i, trip) in enumerate(trips)
        label = trip[1]*trip[2]
        if !(trip in outages)
            draw(
                g;
                outages=outages,
                trip=trip,
                fig=fig[1, 1][(i-1)÷splitlength, (i-1)%splitlength],
                node_labels=nothing,
                edge_labs=nothing,
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
                edge_labs=nothing,
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
        h = deepcopy(g)
        _openings = [i for (i, v_branch) in enumerate(value.(model[:v_branch])) if v_branch == 1]
        _edges = collect(edges(h))
        for (i, f) in enumerate(value.(model[:c_flows]))
            br = e_index_for(h, _edges[i])
            br.p = f / 100
        end
        for v in vertices(h)
            h[label_for(g, v)] = 1 / 100
        end
        h[orig_bus] = -(nv(h) - 1) / 100

        draw(h, outages=openings, layout=layout, title="Connectivity flows")
    end
end

function warmstart!(tomodel, frommodel, var)
    set_start_value.(tomodel[var], value.(frommodel[var]))
end

function retrieve_flows!(g::MetaGraph, model::Model, case)
    !is_solved_and_feasible(model) && @error "model not solved"
    for br in eachindex(model[:flows][case, :, :])
        g[br...].p = value(model[:flows][case, br...])
    end
end

function openings(model::Model)
    is_solved_and_feasible(model) ?
    [index for index in eachindex(model[:v_branch]) if value(model[:v_branch][index]) == 1] :
    ELabel[]
end

function draw_TNR_result(g, model, r; coord=nothing, slack_buses=nothing, slack_phases=nothing, eq::Union{Equivalent,Nothing}=nothing)
    fig = Figure(size=(1500, 800), fontsize=20)

    g_base=copy(g)
    if isnothing(slack_buses)
        dc_flow_optim!(g_base, equivalent=eq)
    else
        dc_flow_optim!(g_base, 
            slack_buses=slack_buses,
            slack_phases=slack_phases, equivalent=eq)
    end
    layout = isa(coord, Dict{String,Point2}) ? coord : Stress(Ptype=Float32)
    draw(g_base, fig=fig[1, 1], layout=layout, title="$(r.tnr_pf)")

    _openings = ELabel[]
    applied_conf = BusConf[]

    display(solution_summary(model))

    if is_solved_and_feasible(model)
        store_result(model, "store_result.csv")
        value(model[:overload]) > 0 && @info "Overload: $(value(model[:overload]))"

        if r.allow_branch_openings || r.OTS_only
            _openings = openings(model)
            println("Open: $_openings")
        end

        if r.bus_confs ≠ BusConf[] && !r.OTS_only
            println("v_bus: $(value.(model[:v_bus]))")
            applied_conf = BusConf[r.bus_confs[i] for (i, v) in enumerate(value.(model[:v_bus])) if isapprox(v, 1, atol=1e-9)]
            println("Applied_conf: $applied_conf")
        end
        # g_result, openings_result = add_subBus(g, applied_conf, _openings) TODO: to be updated with Bus confs !
        g_result = copy(g)
        retrieve_flows!(g_result, model, n_case(r))
        # dc_flow!(g_result, outages=[openings_result])
        draw(g_result, fig=fig[1, 2], outages=_openings, layout=layout, title="objective value:$(round(objective_value(model); digits = 5 ))")
    else
        @warn "NOT FEASIBLE"
    end

    display(fig)
end

function draw_model_case(g::MetaGraph, model::Model, case; kwargs...)
    h = copy(g)
    retrieve_flows!(h, model, case)
    draw(h; kwargs...)
end

function branch_status(model::Model)
    !is_solved_and_feasible(model) && return nothing
    value.(model[:v_branch]) .== 1
end

function display_bus_and_edge_labels(g::MetaGraph)
    println("Bus Labels:")
    for v in 1:nv(g)
        println(v, ":", label_for(g, v))
    end
    println("\nEdge Labels:")
    for (i, elab) in enumerate(edge_labels(g))
        println(i, ":", elab)
    end
end

function display_all_values(model::Model)
    for v in all_variables(model)
        println("$v: $(value(v))")
    end
end

function sandbox(r, equivalent)
    m, r = init_model(r)
    _fixed_buses = Int[r.bus_orig_id]
    _fixed_phases = Float64[0.0]
    create_variables!(m, r, _fixed_buses, equivalent)
    c_w!(m, r)
    # @constraint(m, [c in cases(r)], m[:c_w][c,:] .== 0)

    OTS_flows_phases!(m, r, _fixed_buses, _fixed_phases, equivalent)
    OTS_N_connectedness!(m, r, _fixed_buses, equivalent)

    OTS_N_1_connectedness!(m, r, _fixed_buses, equivalent)

    N_balance!(m, r, _fixed_buses, equivalent)
    OTS_balance!(m, r, _fixed_buses, equivalent)

    loadloss!(m, r, _fixed_buses, equivalent)

    overload!(m, r, equivalent)

    # # branch_status!(m, branch_status)

    @objective(m, Min,
        sum(m[:lostload][c] for c in n_1cases(r)))


    # @constraint(m, m[:v_branch_e][2] == 1)
    # @constraint(m, m[:v_branch][3] == 1)


    optimize!(m)
    if is_solved_and_feasible(m)
        @info "Model OK"
    else
        @warn "Model Infeasible"
    end
    m, r
end

function filter_model(model::Model, st::String)
    for line in eachsplit(string(model), "\n")
        occursin(st, line) && println(line)
    end
end

