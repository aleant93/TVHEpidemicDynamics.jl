"""
    uniform(h::Hypergraph, α::Float64; kwargs...)

 Select a random sample of `α` nodes or `α` hyperdeges of the
 hypergraph `h` to immunize.
"""
function uniform(h::Hypergraph, α::Float64; kwargs...)
    rng = MersenneTwister(1234)
    to_return = shuffle!(rng, collect(1:nhv(h)))
    to_return[1:ceil(Int, length(to_return)*α)]
end



"""
    degrees(h::Hypergraph, α::Float64; kwargs...)

 Select `α` nodes or `α` hyperdeges of the hypergraph `h` to immunize
 according to their degree in `h`.

 The algorithm selects `α` nodes with the higher degree.
"""
function degrees(h::Hypergraph, α::Float64; kwargs...)
    d = [length(gethyperedges(h, v)) for v in 1:nhv(h)]

    to_return = sort!(collect(1:nhv(h)), by = x -> d[x], rev = true)
    to_return[1:ceil(Int, length(to_return)*α)]
end



"""
    centrality(h::Hypergraph, α::Float64; kwargs...)

 Select `α` nodes or `α` hyperdeges of the hypergraph `h` to immunize
 according to their betweeness centrality in `h`.

 The algorithm selects `α` nodes with the higher bc values.
"""
function centrality(h::Hypergraph, α::Float64; kwargs...)
    m = SimpleHypergraphs.adjacency_matrix(h; s=2)

    g = LightGraphs.SimpleGraph(m)
    bc = LightGraphs.betweenness_centrality(g)

    to_return = sort!(collect(1:nhv(h)), by = x -> bc[x], rev = true)
    to_return[1:ceil(Int, length(to_return)*α)]
end



"""
    random_walk(h::Hypergraph, α::Float64; kwargs...)

 Select `α` nodes or `α` hyperdeges of the hypergraph `h` to immunize
 according to their random walk centrality in `h`.

 The algorithm selects `α` nodes with the higher rw values.
"""
function random_walk(h::Hypergraph, α::Float64; kwargs...)
    rwc = Dict{Int, Int}(v => 0 for v in 1:nhv(h))

    for r = 1:100
        for v = 1:nhv(h)
            length(gethyperedges(h, v)) == 0 && continue
            rwc[SimpleHypergraphs.random_walk(h, v)] += 1
        end
    end

    to_return = sort!(collect(1:nhv(h)), by = x -> rwc[x], rev = true)
    to_return[1:ceil(Int, length(to_return)*α)]
end



"""
    acquaintance(h::Hypergraph, α::Float64; kwargs...)

 Select `α` nodes or `α` hyperdeges of the hypergraph `h` to immunize
 according to the acquaintance strategy. It also makes use of local
 knowledge selecting the neighbor with higher degree.
"""
function acquaintance(h::Hypergraph, α::Float64; kwargs...)
    rng = MersenneTwister(1234)
    nodes = shuffle!(rng, collect(1:nhv(h)))

    n_to_immunize = ceil(Int, length(nodes)*α)
    to_immunize = Set{Int}()

    for n in nodes
        neighbors = Set{Int}()
        for he in keys(h.v2he[n])
            union!(neighbors, keys(h.he2v[he]))
        end
        delete!(neighbors, n) #remove v from its neighborhood

        #consider only the nodes that have at least one neighbor
        length(neighbors) == 0 && continue

        # we want to select the neighbor with higher degree
        n_degrees = [(length(gethyperedges(h, n)), n) for n in neighbors]
        chosen = maximum(n_degrees)[2]

        push!(to_immunize, chosen)

        length(to_immunize) == n_to_immunize && break
    end

    return collect(to_immunize)
end



"""
    lockdown(h::Hypergraph, α::Union{Int, Float64}; kwargs...)

 Close all location indicated in `kwargs[:path]`.
"""
function lockdown(h::Hypergraph, α::Union{Int, Float64}; kwargs...)
    return deserialize(kwargs[:path])
end
