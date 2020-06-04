using TVHEpidemicDynamics
using SimpleHypergraphs
using PyPlot

function hasNeighbours(h, node)
    for he in keys(gethyperedges(h, parse(Int, node)))
        if length(getvertices(h, he)) >= 2
            return true
        end
    end
    return false
end

"""
Select a fraction `p` of nodes and immunize one of their neighbours
"""
dataset = "/home/antonio/Scrivania/Tesi/data/dataset_TSMC2014_TKY.txt"
header = [:userid, :venueid, :catid, :catname, :lat, :lng, :timezoneoffset, :UTCtime]
dateformat = "e u d H:M:S +0000 Y" # y-m-dTH:M:SZ

Δₕ = 4
Δₘ = convert(Dates.Millisecond, Dates.Hour(Δₕ))

δₛ = 60
δₘ = convert(Dates.Millisecond, Dates.Minute(δₛ))

df, intervals, node_index_map, he_index_map =
    buildparamsdirect(
        dataset, header, dateformat;
        Δ = Δₘ, δ = δₘ)

h = inithg(df, node_index_map, he_index_map)

nodesWithNeighbours = filter(node->hasNeighbours(h, node), collect(keys(node_index_map)))
pvalues = [0.15, 0.20, 0.25, 0.30]
for p in pvalues
    limit = trunc(Int, length(nodesWithNeighbours) * p)
    immunized = []

    # Search in a node hyperedges
    i = 1
    while length(immunized) <= limit && i <= length(nodesWithNeighbours)
        current = nodesWithNeighbours[i]
        hyperedges = collect(keys(gethyperedges(h, parse(Int, current))))
        j = 1
        toInsert = true
    # a neighbour not already in `immmunized`, if any
        while (j <= length(hyperedges) && toInsert)
            neighbours = collect(keys(getvertices(h, hyperedges[j])))
            filtered = filter(neighbour->!(neighbour in immunized), neighbours)
            if (!isempty(filtered))
                push!(immunized, filtered[1])
                toInsert = false
            end
            j += 1
        end
        i += 1
    end

    infected_distribution = Vector{Array{Float64,1}}()

    linestyles = ["solid", "dashed", "dashdot", "dotted"]
    markers = ["", "", "", "", "x", "+"]

    linestyle = 1
    marker = 1
    clf()
    figure(figsize = (7, 4))
    labels = Array{String,1}()
    mypath = "/home/antonio/Projects/TVHEpidemicDynamics.jl/src/experiments/results/acquaintance/csv/"
    mytitle = "$p"
    plotpath = "/home/antonio/Projects/TVHEpidemicDynamics.jl/src/experiments/results/acquaintance/plot/$(mytitle).png"

    SIS_per_infected_sim =
        TVHSIS(
            df,
            intervals,
            node_index_map,
            he_index_map,
            4,# Δ
            convert(Dates.Millisecond, Dates.Minute(1));# δ
            path = mypath,
            mytitle = mytitle,
            node_imm_strategy = node->(node in immunized),
            he_imm_strategy = he->false
            )

    for t = 1:length(keys(SIS_per_infected_sim)) - 1
        push!(infected_distribution, get(SIS_per_infected_sim, t, Array{Float64,1}()))
    end
    
    ylim(bottom = 0.0, top = 0.6)
    plot(infected_distribution, linestyle = linestyles[linestyle], marker = markers[marker], markevery = 10, markersize = 6.5)
    
    ylabel("Infected nodes in %", fontweight = "semibold", fontsize = "x-large", labelpad = 10)
    
    xlabel("Time intervals", fontweight = "semibold", labelpad = 10, fontsize = "x-large")
    
    tick_params(labelsize = "large")
    
    push!(labels, "acquaintance $p")
    linestyle = (linestyle + 1) % (length(linestyles) + 1)
    marker = (marker + 1) % (length(markers) + 1)
    
    if linestyle == 0
        linestyle = 1
    end
    if marker == 0
        marker = 1
    end
    
    legend(labels, fontsize = "large", ncol = 2)
    plt.tight_layout(.5)
    savefig(plotpath)
end