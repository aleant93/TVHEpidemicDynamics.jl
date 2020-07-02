"""
    uniform(h::Hypergraph)

Select a random set of nodes to immunize
"""
function uniform(h::Hypergraph)
    n = nhv(h)
    Random.seed!(0)
    random_nodes = shuffle(1:n)
    random_nodes
end

"""
    contacts_based(h::Hypergraph)

Select the nodes with the highest contacts number to immunize
"""
function contacts_based(h::Hypergraph)
    # Count nodes in the same hyperedges of `v`
    function countcontacts(h, v)
        contacts = []
        for he in keys(gethyperedges(h, v))
            neighbours = collect(keys(getvertices(h, he)))
            contacts = union(contacts, neighbours)
        end
        length(contacts)
    end

    sorted = sort(1:nhv(h), by=x -> countcontacts(h, x), rev=true)
    sorted
end

"""
    centrality(h::Hypergraph)

Select the nodes with the highest centrality to immunize
"""
function centrality(h)
    centrality_dict = Dict{Int,Int}()
    n = nhv(h)
    map(node -> centrality_dict[node] = 0, 1:n)
    rounds = 10
    for i = 1:rounds
        for node = 1:n
            # A node may not be in no hyperedge
            try
                centrality_dict[random_walk(h, node)] += 1
            catch
                continue
            end
        end
    end
    sorted = sort(collect(centrality_dict), by=x -> x[2], rev=true)
    to_immunize = map(x -> x[1], sorted)
    to_immunize
end

"""
    acquaintance(h::Hypergraph)

Select a fraction of nodes with at least a neighbour,
choose one of their neighbours to immunize
"""
function acquaintance(h::Hypergraph)
    function hasNeighbours(h, node)
        for he in keys(gethyperedges(h, node))
            if length(getvertices(h, he)) >= 2
                return true
            end
        end
        return false
    end

    # select nodes with at least a neighbour
    nodesWithNeighbours = filter(node -> hasNeighbours(h, node), 1:nhv(h))

    to_immunize = []
    Random.seed!(0)
    for node in nodesWithNeighbours
        hyperedges = collect(keys(gethyperedges(h, node)))
        neighbours = []
        for he in hyperedges
            neighbours = union(neighbours, collect(keys(getvertices(h, he))))
        end
        filtered = filter(neighbour -> !(neighbour in to_immunize), neighbours)
        if (!isempty(filtered))
            push!(to_immunize, filtered[rand(1:length(filtered))])
        end
    end

    to_immunize
end
