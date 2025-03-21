using DrWatson
@quickactivate("tnr")

include("../src/GraphUtils.jl")

@proto struct OTS_cc
    buses::Vector{Int}
    l::Float64
    g::Float64
end

function Base.show(io::IO, occ::OTS_cc)
    print("OTS_cc: Buses:", occ.buses, "  load:", occ.l, "  gen:", occ.g)
end
@proto mutable struct Equivalent
    bus_labels::Vector{String}
    buses::Vector{Int} # ids of the connecting bus when connecting to a target network
    p::Vector{Float64}
    gen_slope::Vector{Float64}
    B::Matrix{Float64}
    cc::Vector{OTS_cc}
    cc_id::Vector{Int}
    allow_connector_opening::Bool
    f::Vector{Float64}
    f_max::Vector{Float64}
    F::Matrix{Float64}
end

function Equivalent(
    bus_labels::Vector{String},
    buses::Vector{Int},
    p::Vector{Float64},
    gen_slope::Vector{Float64},
    B::Matrix{Float64},
    cc::Vector{OTS_cc},
    f::Vector{Float64}=Float64[],
    f_max::Vector{Float64}=Float64[],
    F::Matrix{Float64}=Float64[;;];
    allow_connector_opening::Bool=false,
)

    cc_id = zeros(Int, length(buses))
    for (id, _cc) in enumerate(cc)
        foreach(bus_id -> cc_id[bus_id] = id, indexin(_cc.buses, buses))
    end
    Equivalent(bus_labels, buses, p, gen_slope, B, cc, cc_id, allow_connector_opening, f, f_max, F)
end

function Base.show(io::IO, e::Equivalent)
    println(io, "Equivalent: buses: ", e.bus_labels, " -> ", e.buses)
    println(io, "  p: ", e.p)
    println(io, "  gen_slope: ", e.gen_slope)
    println(io, "  B: ", e.B)
    println(io, "  allow_connector_opening: ", e.allow_connector_opening)
    println(io, "  cc:")
    foreach(cc -> println(io, "    ", cc), e.cc)
    if !isempty(e.f)
        println("  flows")
        println(io, "    f:", e.f)
        println(io, "    f_max:", e.f_max)
        println(io, "    F:", e.F)
    end
end

get_bus_id(eq::Equivalent, bus::Int) = findfirst(==(bus), eq.buses)

get_cci(eq::Equivalent, bus::Int) = eq.cc_id[get_bus_id(eq, bus)]