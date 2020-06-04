using TVHEpidemicDynamics
using PyPlot

"""
Select a random set of nodes and vaccinate them
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

percentages = [0.70, 0.75, 0.80, 0.85, 0.90]
sorted = sort(collect(node_index_map), by = x->x[2])
# sorted = sort(collect(he_index_map), by = x->x[1])

# n% of nodes/hyperedges are infectious
for n in percentages

    index = trunc(Int, length(sorted) * n)
    infected_nodes = map(x->parse(Int64, x[1]), sorted[1:index])
    # infected_hes = map(x->x[2], sorted[1:index])

    infected_distribution = Vector{Array{Float64,1}}()

    linestyles = ["solid", "dashed", "dashdot", "dotted"]
    markers = ["", "", "", "", "x", "+"]

    linestyle = 1
    marker = 1
    clf()
    figure(figsize = (7, 4))
    labels = Array{String,1}()
    mypath = "/home/antonio/Projects/TVHEpidemicDynamics.jl/src/experiments/results/uniform/nodes/csv/"
    mytitle = "$n"
    plotpath = "/home/antonio/Projects/TVHEpidemicDynamics.jl/src/experiments/results/uniform/nodes/plot/$(mytitle).png"
    # plotpath = "/home/antonio/Projects/TVHEpidemicDynamics.jl/src/experiments/results/uniform/hes/plot/$(mytitle).png"

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
            # node_imm_strategy = node->false,
            # he_imm_strategy = he->!(he in infected_hes)
            )

    for t = 1:length(keys(SIS_per_infected_sim)) - 1
        push!(infected_distribution, get(SIS_per_infected_sim, t, Array{Float64,1}()))
    end

    ylim(bottom = 0.0, top = 0.6)
    plot(infected_distribution, linestyle = linestyles[linestyle], marker = markers[marker], markevery = 10, markersize = 6.5)

    ylabel("Infected nodes in %", fontweight = "semibold", fontsize = "x-large", labelpad = 10)

    xlabel("Time intervals", fontweight = "semibold", labelpad = 10, fontsize = "x-large")

    tick_params(labelsize = "large")

    push!(labels, "uniform $n")
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