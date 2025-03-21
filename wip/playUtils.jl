include("../src/DrawGrid.jl")
include("../src/GraphUtils.jl")
include("../src/dcpf.jl")

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
        add_constraint(g, b -> b.p_max = 0.25)
        g["12", "13"].p_max = 0.4
        g["6", "8"].p_max = 0.4
        g["4", "6"].p_max = 0.4
    end


    balance!(g, all_non_zero_uniform)

    g, bus_confs, coord
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

function calcanddraw(g; bus_orig=nothing,
    outages::Union{AbstractArray{Int},AbstractArray{Tuple{String,String}}}=Int[], trip::Union{Nothing,Int,Tuple{String,String}}=nothing,
    label="", eq=nothing, kwargs...)
    _outages = outages isa AbstractArray{Tuple{String,String}} ? e_code_for(g, outages) : outages

    bus_orig = isnothing(bus_orig) ? label_for(g, 1) : bus_orig
    _trip = trip isa Tuple{String,String} ? e_code_for(g, trip...) : trip
    g_cc = connected(g, bus_orig; outages=_outages, trip=_trip, eq=eq)

    if nv(g) == nv(g_cc)
        h = g
        houtages = _outages
        htrip = _trip
    else
        h = g_cc
        balance!(h)
        houtages = Int[]
        htrip = nothing
        label = label == "" ? "" : L"\textbf{%$label}(%$(round(Int, 100*(total_load(g)-total_load(g_cc)))))"
    end
    dc_draw = copy(h)
    dc_flow_optim!(dc_draw; outages=houtages, trip=htrip, equivalent=eq)
    draw(dc_draw, outages=houtages, trip=htrip, title=label,
        edge_width=br -> br.outage || br.trip ? 6 : 2,
        ; kwargs...)
end

function griddraw(g, bus_orig::String, trips::Union{AbstractArray{Int},AbstractArray{Tuple{String,String}}}, outages::Union{AbstractArray{Int},AbstractArray{Tuple{String,String}}}=Int[], eq=nothing; kwargs...)
    _trips = isa(trips, AbstractArray{Tuple{String,String}}) ? e_code_for(g, trips) : trips
    _outages = isa(outages, AbstractArray{Tuple{String,String}}) ? e_code_for(g, outages) : outages
    fig = Figure(size=(600, 500), fontsize=20)
    splitlength = 5
    all_edges = collect(edges(g))
    for (i, trip) in enumerate(_trips)
        e = all_edges[trip]
        label = "$(label_for(g, src(e)))$(label_for(g, dst(e)))"

        if !(trip in _outages)
            draw(
                g;
                outages=_outages,
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
                outages=_outages,
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

function retrieve_flows!(g::MetaGraph, model::Model, case::Int)
    !is_solved_and_feasible(model) && @error "model not solved"
    for (i, flow) in enumerate(value.(model[:flows][case, :]))
        g[e_label_for(g, i)...].p = flow
    end
end

function openings(model::Model)
    is_solved_and_feasible(model) ?
    [i for (i, v_branch) in enumerate(value.(model[:v_branch])) if v_branch == 1] :
    Int[]
end

function draw_TNR_result(g, model, r; coord=nothing, fixed_buses=nothing, fixed_phases=nothing, eq::Union{Equivalent,Nothing}=nothing)
    fig = Figure(size=(1500, 800), fontsize=20)

    if isnothing(fixed_buses)
        g_base = dc_flow(g, equivalent=eq)
    else
        g_base = dc_flow(g;
            fixed_buses=isa(fixed_buses, AbstractArray{Int}) ? fixed_buses : to_code(g, fixed_buses),
            fixed_phases=fixed_phases, equivalent=eq)
    end
    layout = isa(coord, Dict{String,Point2}) ? coord : Stress(Ptype=Float32)
    draw(g_base, fig=fig[1, 1], layout=layout, title="$(r.tnr_pf)")

    _openings = Int[]
    applied_conf = BusConf[]

    display(solution_summary(model))

    if is_solved_and_feasible(model)
        store_result(model, "store_result.csv")
        value(model[:overload]) > 0 && @info "Overload: $(value(model[:overload]))"

        if r.allow_branch_openings || r.OTS_only
            _openings = openings(model)
            println("Open: $_openings")
            println([collect(edge_labels(g))[i] for i in _openings])
        end

        if r.bus_confs ≠ BusConf[] && !r.OTS_only
            println("v_bus: $(value.(model[:v_bus]))")
            applied_conf = BusConf[r.bus_confs[i] for (i, v) in enumerate(value.(model[:v_bus])) if isapprox(v, 1, atol=1e-9)]
            println("Applied_conf: $applied_conf")
        end
        g_result, openings_result = add_subBus(g, applied_conf, _openings)
        retrieve_flows!(g_result, model, n_case(r))
        # dc_flow!(g_result, outages=[openings_result])
        draw(g_result, fig=fig[1, 2], outages=openings_result, layout=layout, title="objective value:$(round(objective_value(model); digits = 5 ))")
    else
        @warn "NOT FEASIBLE"
    end

    display(fig)
end

function draw_model_case(g::MetaGraph, model::Model, case::Int; kwargs...)
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
    for v in 1:ne(g)
        println(v, ":", e_label_for(g, v))
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