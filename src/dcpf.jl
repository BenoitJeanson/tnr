using SparseArrays
using MetaGraphsNext

function dc_flow!(g::MetaGraph; outages = [], trip = nothing)
    A = incidence_matrix(g; oriented = true)'
    for outage in outages
        A[outage,:] .= 0
    end
    if trip ≠ nothing
        A[trip,:] .= 0
    end
    D = spdiagm([e_index_for(g,e).b for e in edges(g)])
    B = A'*D*A
    B_inv = inv(Matrix(B + fill(1/nv(g), nv(g), nv(g)))) - fill(1/nv(g), nv(g), nv(g))
    
    p_to_f = D*A*B_inv
    flows = p_to_f * [g[label_for(g, v)] for v in vertices(g)]
    for (i,e) in enumerate(edges(g))
        e_index_for(g, e).p = flows[i]
    end
    g
end

function dc_flow(g::MetaGraph; outages = [], trip = nothing)
    h=deepcopy(g)
    dc_flow!(h, outages = outages, trip = trip)
end
