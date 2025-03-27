include("../src/GraphUtils.jl")

@proto struct OTS_cc
    buses::Vector{VLabel}
    l::Float64
    g::Float64
end

function Base.show(io::IO, occ::OTS_cc)
    print("OTS_cc: Buses:", occ.buses, "  load:", occ.l, "  gen:", occ.g)
end
@proto mutable struct Equivalent
    buses::Vector{VLabel}
    p::Vector{Float64}
    gen_slope::Vector{Float64}
    B::Matrix{Float64}
    cc::Vector{OTS_cc}
    cc_id::Vector{Int}
    allow_connector_opening::Bool
    internal_branches::Vector{ELabel}
    f::Vector{Float64}
    f_max::Vector{Float64}
    F::Matrix{Float64}
end

function Equivalent(
    buses::Vector{VLabel},
    p::Vector{Float64},
    gen_slope::Vector{Float64},
    B::Matrix{Float64},
    cc::Vector{OTS_cc},
    internal_branches::Vector{ELabel}=ELabel[],
    f::Vector{Float64}=Float64[],
    f_max::Vector{Float64}=Float64[],
    F::Matrix{Float64}=Float64[;;];
    allow_connector_opening::Bool=false,
)

    cc_id = zeros(Int, length(buses))
    for (id, _cc) in enumerate(cc)
        foreach(bus_id -> cc_id[bus_id] = id, indexin(_cc.buses, buses))
    end
    Equivalent(buses, p, gen_slope, B, cc, cc_id, allow_connector_opening, internal_branches, f, f_max, F)
end

function Base.show(io::IO, e::Equivalent)
    println(io, "Equivalent: buses: ", e.buses)
    println(io, "  p: ", e.p)
    println(io, "  gen_slope: ", e.gen_slope)
    println(io, "  B: ", e.B)
    println(io, "  allow_connector_opening: ", e.allow_connector_opening)
    println(io, "  cc:")
    foreach(cc -> println(io, "    ", cc), e.cc)
    if !isempty(e.f)
        println("  flows")
        println(io, "    internal_branches: ", e.internal_branches)
        println(io, "    f:", e.f)
        println(io, "    f_max:", e.f_max)
        println(io, "    F:", e.F)
    end
end

get_bus_id(eq::Equivalent, bus::VLabel) = findfirst(==(bus), eq.buses)

get_cci(eq::Equivalent, bus::VLabel) = eq.cc_id[get_bus_id(eq, bus)]