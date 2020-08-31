#########################
# Plotting new infected
# per location per day
#
# Produces a heat map representing most infeced locations over time
########################
using Pkg
Pkg.activate(".")
using PyPlot
using Serialization
using TVHEpidemicDynamics

project_path = dirname(pathof(TVHEpidemicDynamics))
output_path = joinpath(project_path, "experiments", "spreading", "ble", "heatmap")

# data_to_plot = infectedXloc_data["exp415-low"]
# sim_data = intervals_data["415"]
#
# serialize("ble_data_heatmap_infectedPerLocation.data", infectedXloc_data)
# serialize("sim_data_infectedPerLocation.data", intervals_data)

data_to_plot = deserialize(joinpath(output_path, "ble_data_heatmap_infectedPerLocation.data"))["exp415-low"]
sim_data = deserialize(joinpath(output_path, "sim_data_infectedPerLocation.data"))["415"]

he2loc = Dict{Int, String}()

for loc in keys(sim_data[:loc2he])
    push!(
        he2loc,
        sim_data[:loc2he][loc] => (loc)
    )
end

#00 01 02 30 31 25 29 08 13 28 24 23 27 26 22 14 15 16 17 18 19 20 21 03 04 06 07 08 09 10 11 12
order = ["00" "01" "02" "30" "31" "25" "29" "08" "13" "28" "24" "23" "27" "26" "22" "14" "15" "16" "17" "18" "19" "20" "21" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12"]
_ids = ["rpi-$(id)" for id in order]

ids = reverse(_ids, dims=2)

_labels = ["00", "01", "02", "30", "31", "25", "29", "08", "13", "28", "24", "23", "27", "26", "22", "14", "15", "16", "17", "18", "19", "20", "21", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]
labels = ["rpi-$(id)" for id in reverse(_labels)]

#days = ["$(index)" for (index,id) in sort!(collect(keys(data_to_plot[1][1])))]

days = sort!(collect(keys(data_to_plot[1][1])))

m = fill(0.0, length(ids), length(days))


for he = 1:length(he2loc)
    for day in days
        val = 0.0
        for iter in 1:10
            val += data_to_plot[iter][he][day]
        end
        val != 0 && (val \= 10) > 0.0

        print(findfirst(x -> x == day, days))
        print( " idLoc ", he2loc[he])
        println( "   indexLoc ", findfirst(x -> x == he2loc[he], ids)[2])

        m[findfirst(x -> x == he2loc[he], ids)[2], findfirst(x -> x == day, days)] = val
    end
end

m

clf()
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)

#im = ax.imshow(m)
im = ax.imshow(m, interpolation="gaussian")

# # We want to show all ticks...
yticks(0:length(ids)-1, labels)
xticks(0:length(days)-1, days, rotation=80)

gcf()

plt.tight_layout(.5)
savefig(joinpath(output_path, "heatmap.png"))
