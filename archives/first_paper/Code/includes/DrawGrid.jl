# module DrawGrid
using LaTeXStrings
using Graphs
using MetaGraphsNext
using CSV
using CairoMakie
using GraphMakie
using GraphMakie.NetworkLayout

include("GraphUtils.jl")

export generate_csv, load_coord, complete_subBus_coords, draw

CairoMakie.activate!(type = "svg")

"""
    generate_csv(Dict{String, Point2}) -> String
    
Generates a String that represents a CSV file from the dict of Points
"""
function generate_csv(dic)
    println("bus_i,x,y")
    foreach(l->println("$l, $(dic[l][1]), $(dic[l][2])"), sort(collect(keys(dic)), by=x-> parse(Int, x)))
end

"""
    load_coord(filename) -> Dict{String, Point2}

Load a csv file containing the coordinates
"""
function load_coord(file, buslabels = num->"$num")
    csv_data = CSV.File(file)
    coord = Dict{String, Point2}()
    foreach(row -> coord[buslabels(row.bus_i)] = Point(row.x, row.y), csv_data)
    coord
end

function complete_subBus_coords(g, init_points)
    res_points = Point2[]
    for lab in labels(g)
        if haskey(init_points, lab)
            push!(res_points, init_points[lab])
        else
            neighbs = [split(label_for(g, c), "-")[1] for c in all_neighbors(g, code_for(g, lab))]
            push!(neighbs, split(lab, "-")[1])
            println(neighbs)
            println(init_points)

            push!(res_points, 1/length(neighbs) * reduce(+, [init_points[lab] for lab in neighbs]))
        end
    end
    res_points
end

function draw(g_orig::MetaGraph;
                fig::Union{Figure, GridPosition, GridSubposition, Nothing} = nothing,
                figposition = Tuple{Int, Int},

                outages = [],
                trip = nothing,
    
                fig_size = (900, 600),
                font_size = 15,
                layout = Stress(Ptype=Float32),
                margin_ratio =7,
                digits = 2,
    
                node_labels = (g, label; digits = 0) -> LaTeXString("\\textbf{$label:}$(round(Int, 100*g[label]))"),
                marker = (g, label) -> g[label]>=0 ? :rect : :circle,
                node_size = (g, label, max_p) -> sqrt(abs(g[label]))/sqrt(max_p)*50,

                edge_labels = (br; digits = 0) -> L"%$(round(Int,100*br.p))",
                edge_coloring = (br; tol=1e-6) -> br.outage || br.trip ? :black : (abs(br.p>br.p_max + tol) ? :red : :green3),
                edge_style = branch -> branch.trip ? (:dot, :dense) : :solid,
                edge_width = br -> br.outage || br.trip ? 4 : 2,
                arrow_size = 25,
                title = ""
                )

    g = copy(g_orig)
    max_p = maximum(collect(labels(g)) .|> l->abs(g[l]))
    gLabels = [label_for(g, v) for v in vertices(g)]
    to_reverse = []
    for (i,e) in enumerate(edges(g))
        br = e_index_for(g, e)
        if i in outages
            br.outage = true
        end
        if i == trip
            br.trip = true
        end
        if br.p<0
            push!(to_reverse, e)
        end
    end
    for e in to_reverse
        b = e_index_for(g, e)
        b.p = -b.p
        rem_edge!(g, src(e), dst(e))
        g[label_for(g,dst(e)), label_for(g,src(e))] = b
    end

    if isnothing(fig)
        fig = Figure(size = fig_size, fontsize = font_size)
        ax = Axis(fig[1, 1], title = title)
    else
        ax = Axis(fig, title = title)
    end

    x_min, x_max, y_min, y_max = reduce(((x_min, x_max, y_min, y_max), pt) 
        -> (min(pt[1], x_min), max(pt[1], x_max), min(pt[2], y_min), max(pt[2], y_max)), 
        layout isa Dict{String, Point2} ? values(layout) : layout(g) , init = (0,0,0,0))
    delta_x, delta_y = ((x_max - x_min) , (y_max - y_min)) ./ margin_ratio
    xlims!(ax, x_min-delta_x, x_max+delta_x)
    ylims!(ax, y_min-delta_y, y_max+delta_y)
    
    g_labels = [label_for(g, v) for v in vertices(g)]
    branches = [e_index_for(g, e) for e in edges(g)]

    graphplot!(ax, g;
        layout = layout isa Dict{String, Point2} ? complete_subBus_coords(g, layout) : layout,
        node_size = isnothing(node_size) ? 0 : [node_size(g, l, max_p) for l in gLabels],
        node_attr = (marker = [marker(g, l) for l in g_labels], color= :white, strokecolor = :red, strokewidth = 3,),
        nlabels = isnothing(node_labels) ? nothing : [node_labels(g, l, digits=digits) for l in g_labels],
        nlabels_align = (:center,:top),
        nlabels_attr=(rotation=0,),
        nlabels_distance=-40,
        elabels = isnothing(edge_labels) ? nothing : [edge_labels(b, digits=digits) for b in branches],
        edge_plottype = :beziersegments,
		edge_attr = (linestyle = [edge_style(b) for b in branches],),
        edge_color = [edge_coloring(b) for b in branches],
        edge_width = [edge_width(b) for b in branches],
        arrow_size = [b.outage || b.trip ? 0 : arrow_size for b in branches],
        arrow_shift = :end,
        )
    hidedecorations!(ax)
    hidespines!(ax)

    ax.aspect = DataAspect()
    fig
end

# end