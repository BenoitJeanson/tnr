using ProtoStructs
using MetaGraphsNext

include("GraphUtils.jl")

@enum TNR_PF_TYPE tnr_pf_pst tnr_pf_phase
@proto struct TNR
    g::MetaGraph
    n_1_connectedness::Bool
    bus_confs::Vector{BusConf}
    allow_branch_openings::Bool
    OTS_only::Bool
    tnr_pf::TNR_PF_TYPE
    opf::Bool
    bus_orig::String

    bus_orig_id::Int
    outages::AbstractArray{Int}
    A::AbstractMatrix{Int}
    branches::Vector{Branch}
    p::Vector{Float64}
    bus_to_conf_ids::Dict{Int,Vector{Int}}
    branch_to_conf_ids::Dict{Int,Vector{Int}}
    bus_branch_to_conf_ids::Dict{Tuple{Int,Int},Vector{Int}}
end

function TNR(g::MetaGraph,
    contingencies::AbstractArray{Int},
    n_1_connectedness::Bool,
    bus_confs::Vector{BusConf},
    allow_branch_openings::Bool,
    OTS_only,
    tnr_pf,
    opf,
    bus_orig)

    outages = collect(contingencies)
    push!(outages, 0)
    A = incidence_matrix(g; oriented=true)'

    TNR(
        g=g,
        n_1_connectedness=n_1_connectedness,
        bus_confs=bus_confs,
        allow_branch_openings=allow_branch_openings,
        OTS_only=OTS_only,
        tnr_pf=tnr_pf,
        opf=opf,
        bus_orig=bus_orig, bus_orig_id=code_for(g, bus_orig),
        outages=outages,
        A=A, branches=[e_index_for(g, e) for e in edges(g)],
        p=[g[label_for(g, v)] for v in vertices(g)], bus_to_conf_ids=reduce((acc, (index, conf)) -> (
                push!(get!(acc, conf.bus, Vector{Int}()), index);
                acc
            ), enumerate(bus_confs),
            init=Dict{Int,Vector{Int}}()
        ), branch_to_conf_ids=reduce((acc, (index, conf)) -> (
                foreach(br -> push!(get!(acc, br, Vector{Int}()), index),
                    (br for sb in conf.subBuses for br in sb.branch_ids));
                acc
            ), enumerate(bus_confs),
            init=Dict(i => Int[] for i in 1:ne(g))
        ), bus_branch_to_conf_ids=reduce(
            (acc, (bc_id, bus, branch)) -> (
                push!(get!(acc, (bus, branch), Vector{Int}()), bc_id);
                acc),
            ((bc_id, bc.bus, branch)
             for (bc_id, bc) in enumerate(bus_confs)
             for subBus in bus_confs[bc_id].subBuses
             for branch in subBus.branch_ids),
            init=Dict{Tuple{Int,Int},Vector{Int}}())
    )
end

nb_cases(tnr::TNR) = length(tnr.outages)

cases(tnr::TNR) = cases = 1:nb_cases(tnr)

n_1cases(tnr::TNR) = 1:nb_cases(tnr)-1

n_case(tnr::TNR) = nb_cases(tnr)

nb_buses(tnr::TNR) = nv(tnr.g)

buses(tnr::TNR) = 1:nb_buses(tnr)

edge_ids(tnr::TNR) = 1:ne(tnr.g)

nb_bus_confs(tnr::TNR) = length(tnr.bus_confs)

bus_conf_ids(tnr::TNR) = 1:nb_bus_confs(tnr)

branch_side2bus(tnr::TNR, br::Int, side::Int) = findfirst(i -> tnr.A[br, i] == (side == 1 ? -1 : 1), 1:nv(tnr.g))

branch_bus2side(tnr::TNR, br::Int, v::Int) = tnr.A[br, v] == 1 ? 2 : -tnr.A[br, v]

is_outage(tnr::TNR, i::Int, j::Int) = j == tnr.outages[i] ? 1 : 0

bus_nb_confs(tnr::TNR, bus) = haskey(tnr.bus_to_conf_ids, bus) ? length(tnr.bus_to_conf_ids[bus]) : 0

incident(tnr::TNR, v::Int) = [e for e in edge_ids(tnr) if tnr.A[e, v] ≠ 0]

opposite(tnr::TNR, e::Int, v::Int) = findfirst(i -> (i ≠ v && tnr.A[e, i] ≠ 0), buses(tnr))

function Base.show(io::IO, t::TNR)
    println(io, "TNR: Graph with ", nv(t.g),
        " buses and ", ne(t.g), " edges, ",
        length(t.outages), " cases. ",
        t.OTS_only ? "Branches opening" :
        (t.allow_branch_openings ? "Buses splitting and branches opening" : "Buses splitting"),".")
end