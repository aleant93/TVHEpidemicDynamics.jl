"""
Select a random set of nodes to immunize
"""
function uniform(h, perc)
    n = nhv(h)
    Random.seed!(0)
    random_nodes = shuffle(1:n)
    limit = trunc(Int, n * perc)
    random_nodes[1:limit]
end

"""
Select the nodes with the highest contacts number to immunize
"""
function contacts_based(h, perc)
    # Count nodes in the same hyperedges of `v`
    function countcontacts(h, v)
        contacts = []
        for he in keys(gethyperedges(h, v))
            neighbours = collect(keys(getvertices(h, he)))
            contacts = union(contacts, neighbours)
        end
        length(contacts)
    end

    n = nhv(h)
    sorted = sort(1:n, by=x->countcontacts(h, x), rev=true)
    limit = trunc(Int, n * perc)
    sorted[1:limit]
end

"""
Select the nodes with the highest centrality to immunize
"""
function centrality(h, perc)
    centrality_dict = Dict{Int,Int}()
    n = nhv(h)
    map(node->centrality_dict[node] = 0, 1:n)
    # for i = 1:n
    #     map(node->centrality_dict[random_walk(h, node)] += 1, 1:nhv(h))
    # end
    for i = 1:n
        for node = 1:n
            try
                centrality_dict[random_walk(h, node)] += 1
            catch
                continue
            end
        end
    end
    sorted = sort(collect(centrality_dict), by=x->x[2], rev=true)
    print(sorted)
    limit = trunc(Int, n * perc)
    to_immunize = map(x->x[1], sorted[1:limit])
    to_immunize
end

"""
Select a fraction of nodes with at least a neighbour,
choose one of their neighbours to immunize
"""
function acquaintance(h, perc)
    function hasNeighbours(h, node)
        for he in keys(gethyperedges(h, node))
            if length(getvertices(h, he)) >= 2
                return true
            end
        end
        return false
    end

    nodesWithNeighbours = filter(node->hasNeighbours(h, node), 1:nhv(h))
    limit = trunc(Int, length(nodesWithNeighbours) * perc)
    if limit == 0
        return []
    end

    # Search in a node hyperedges
    to_immunize = []
    i = 1
    Random.seed!(0)
    while length(to_immunize) <= limit && i <= length(nodesWithNeighbours)
        current = nodesWithNeighbours[i]
        hyperedges = collect(keys(gethyperedges(h, current)))
        j = 1
        toInsert = true
        # a neighbour not already in `to_immmunize`, if any
        while (j <= length(hyperedges) && toInsert)
            neighbours = collect(keys(getvertices(h, hyperedges[j])))
            random_neighbours = shuffle(collect(neighbours))
            filtered = filter(neighbour->!(neighbour in to_immunize), random_neighbours)
            if (!isempty(filtered))
                push!(to_immunize, filtered[1])
                toInsert = false
            end
            j += 1
        end
        i += 1
    end

    to_immunize
end
