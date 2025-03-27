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
    bus_orig::VLabel

    contingencies::AbstractArray{ELabel}
    n_case::String
    n_1cases::Vector{String}
    cases::Vector{String}
    # A::AbstractMatrix{Int}
    # bus_to_conf_ids::Dict{Int,Vector{Int}}
    # branch_to_conf_ids::Dict{Int,Vector{Int}}
    # bus_branch_to_conf_ids::Dict{Tuple{Int,Int},Vector{Int}}
end

branch_to_case(case::ELabel) = case[1] * case[2]

function TNR(g::MetaGraph,
    contingencies::AbstractArray{ELabel},
    n_1_connectedness::Bool,
    bus_confs::Vector{BusConf},
    allow_branch_openings::Bool,
    OTS_only,
    tnr_pf,
    opf,
    bus_orig)

    A = incidence_matrix(g; oriented=true)'
    n_case = "-"
    n_1cases=collect(map(branch_to_case, contingencies))

    TNR(
        g=g,
        n_1_connectedness=n_1_connectedness,
        bus_confs=bus_confs,
        allow_branch_openings=allow_branch_openings,
        OTS_only=OTS_only,
        tnr_pf=tnr_pf,
        opf=opf,
        bus_orig=bus_orig,
        n_case=n_case,
        contingencies=contingencies,
        n_1cases=n_1cases,
        cases=[n_1cases; n_case],
        # A=A,
        # bus_to_conf_ids=reduce((acc, (index, conf)) -> (
        #         push!(get!(acc, conf.bus, Vector{Int}()), index);
        #         acc
        #     ), enumerate(bus_confs),
        #     init=Dict{Int,Vector{Int}}()
        # ),
        # branch_to_conf_ids=reduce((acc, (index, conf)) -> (
        #         foreach(br -> push!(get!(acc, br, Vector{Int}()), index),
        #             (br for sb in conf.subBuses for br in sb.branch_ids));
        #         acc
        #     ), enumerate(bus_confs),
        #     init=Dict(i => Int[] for i in 1:ne(g))
        # ),
        # bus_branch_to_conf_ids=reduce(
        #     (acc, (bc_id, bus, branch)) -> (
        #         push!(get!(acc, (bus, branch), Vector{Int}()), bc_id);
        #         acc),
        #     ((bc_id, bc.bus, branch)
        #      for (bc_id, bc) in enumerate(bus_confs)
        #      for subBus in bus_confs[bc_id].subBuses
        #      for branch in subBus.branch_ids),
        #     init=Dict{Tuple{Int,Int},Vector{Int}}())
    )
end

TNR(g::MetaGraph, bus_orig) = TNR(g, collect(edge_labels(g)), false, BusConf[], true, true, tnr_pf_phase, false, bus_orig)

n_case(tnr::TNR) = tnr.n_case

branch_to_case(tnr::TNR, case::ELabel) = case[1] * case[2]

cases(tnr::TNR) = tnr.cases

n_1cases(tnr::TNR) = tnr.n_1cases

nb_buses(tnr::TNR) = nv(tnr.g)

buses(tnr::TNR) = labels(tnr.g)

branches(tnr::TNR) = edge_labels(tnr.g)

nb_bus_confs(tnr::TNR) = length(tnr.bus_confs)

bus_conf_ids(tnr::TNR) = 1:nb_bus_confs(tnr)

branch_side2bus(tnr::TNR, br::Int, side::Int) = findfirst(i -> tnr.A[br, i] == (side == 1 ? -1 : 1), 1:nv(tnr.g))

branch_bus2side(tnr::TNR, br::Int, v::Int) = tnr.A[br, v] == 1 ? 2 : -tnr.A[br, v]

is_outage(tnr::TNR, case::Int, j::Int) = j == tnr.contingencies[case] ? 1 : 0

is_outage(case::String, e::ELabel) = case == branch_to_case(e)

bus_nb_confs(tnr::TNR, bus) = haskey(tnr.bus_to_conf_ids, bus) ? length(tnr.bus_to_conf_ids[bus]) : 0

incident_signed(tnr::TNR, bus) = incident_signed(tnr.g, bus)

incident(tnr::TNR, bus::VLabel) = incident(tnr.g, bus)

opposite(tnr::TNR, e::Int, v::Int) = findfirst(i -> (i ≠ v && tnr.A[e, i] ≠ 0), buses(tnr))

function Base.show(io::IO, tnr::TNR)
    println(io, "TNR: Graph with ", nv(tnr.g),
        " buses and ", ne(tnr.g), " edges, ",
        length(tnr.contingencies), " cases. ",
        tnr.OTS_only ? "Branches opening" :
        (tnr.allow_branch_openings ? "Buses splitting and branches opening" : "Buses splitting"), ".")
end