using TVHEpidemicDynamics
using SimpleHypergraphs
using PyPlot

"""
Random-walk centrality immunization policy
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

# Evaluate centrality for each node
centrality_dict = Dict()
map(node->centrality_dict[node] = 0, collect(keys(node_index_map)))
n = 10
for i in 1:n
    map(node->centrality_dict[string(random_walk(h, parse(Int, node)))] += 1, collect(keys(node_index_map)))
end

# (1-n) percentage of nodes with least centraliy are infectious
sorted = sort(collect(centrality_dict), by = x->x[2])
percentages = [0.10, 0.15, 0.20, 0.25, 0.30]

for n in percentages
    index = trunc(Int, length(sorted) * n)
    infected_nodes = map(x->x[1], sorted[1:index])

    infected_distribution = Vector{Array{Float64,1}}()

    linestyles = ["solid", "dashed", "dashdot", "dotted"]
    markers = ["", "", "", "", "x", "+"]

    linestyle = 1
    marker = 1
    clf()
    figure(figsize = (7, 4))
    labels = Array{String,1}()
    mypath = "/home/antonio/Projects/TVHEpidemicDynamics.jl/src/experiments/results/centrality/csv/"
    mytitle = "$n"
    plotpath = "/home/antonio/Projects/TVHEpidemicDynamics.jl/src/experiments/results/centrality/plot/$(mytitle).png"

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
            node_imm_strategy = node->!(node in infected_nodes),
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

    push!(labels, "centrality $n")
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