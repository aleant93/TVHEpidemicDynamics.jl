using TVHEpidemicDynamics
using DataFrames
using Dates
using CSV
using PyPlot
using Statistics

"""
    Experiments from sections 5.3 and 5.4 of the paper
    *A Design-Methodology for Epidemic Dynamics via Time-Varying Hypergraphs*
        - 5.3 Direct vs indirect contagions
        - 5.4 Modeling the effect of time
"""

#########################
# Loading simulation params
########################

# Simulation parameters
# Section 5.3 - Direct vs Indirect contagions
fparams = "src/experiments/AAMAS20/configs/53/blebluetooth.csv"
output_path = "src/experiments/AAMAS20/results/53"

# Section 5.4 - Modeling the effect of time
#fparams = "src/experiments/AAMAS20/configs/54/aamas54.csv"
#output_path = "src/experiments/AAMAS20/results/54"

paramsdf = CSV.read(
            fparams;
            copycols = true,
            header = [:exp_id, :data, :label, :per_infected, :c,  :Δ, :δ ,:βd, :βᵢ, :βₑ, :γₑ, :γₐ],
            datarow = 2
        )

# just a trick to group together
# all experiments to show in the same plot
test_data = Dict{String, Array{Any, 1}}()
for params in eachrow(paramsdf)
    push!(
        get!(test_data, params[:exp_id], Array{Any, 1}()),
        params
    )
end


#########################
# Generating model data
########################

# Foursqaure dataset
dataset = "data/blebeacon/BLEBeacon-Dataset-master/Check-In Check-Out Report.csv"#"data/dataset_TSMC2014_TKY.txt"
header = [:entry_id, :userid, :UTCtime, :out_time, :venueid]
dateformat = "Y-m-d H:M:S"

# The simulation will consider only the data
# within this time intervals
firstcheckindate = Dates.DateTime("2016-09-15T00:00:00")
lastcheckindate = Dates.DateTime("2016-10-18T00:00:00")

# The choice of the interval within which
# either an indirect (Δ) or direct (δ) contact
# may occur influences the data the
# simulation is run on.
# For this reason, it is necessary to store
# diffent information according to the
# values of both Δ and δ.
intervals = unique(paramsdf, [:Δ, :δ])[!, [:Δ, :δ]]
intervals_data = Dict{String, Dict{Symbol, Any}}()

for i in eachrow(intervals)
    df, intervals, user2vertex, loc2he =
        generate_model_data(
            dataset,
            header,
            :userid,
            :venueid,
            :UTCtime,
            dateformat;
            Δ = convert(Dates.Millisecond, Dates.Hour(i.Δ)),
            δ = convert(Dates.Millisecond, Dates.Minute(i.δ)),
            maxdate = lastcheckindate,
            mindate = firstcheckindate
        )

    push!(
        get!(intervals_data, "$(i.Δ)$(i.δ)", Dict{Symbol, Any}()),
        :df => df,
        :intervals => intervals,
        :user2vertex => user2vertex,
        :loc2he => loc2he,
    )
end



#CHECKIN DISTRIBUTION OVER THE INTERVALS
for key in keys(intervals_data)
    i2c = Array{Int, 1}()
    _intervals = intervals_data[key][:intervals]

    for t in 1:length(_intervals)
        mindate = get(_intervals, t, 0).first
        maxdate = get(_intervals, t, 0).second

        _df = filter(
            r-> (isnothing(mindate) || r[:UTCtime] >= mindate) &&
                (isnothing(maxdate) || r[:UTCtime] <= maxdate),
                intervals_data[key][:df]
        )

        push!(
            i2c,
            nrow(_df)
        )
    end

    push!(
        get!(intervals_data, key, Array{Int, 1}()),
        :checkins => i2c
    )
end


#########################
# Initialization of infected nodes
########################

# For the reproducibility of the experiments,
# the infected nodes at start as to be the same
per_infected = unique(paramsdf, [:per_infected])[!, [:per_infected]]
per_infected_data = Dict{Int, Array{Int, 1}}()

users = keys(intervals_data[collect(keys(intervals_data))[1]][:user2vertex])

for p in eachrow(per_infected)
    vstatus = rand(1:1, 1, length(users)) #size(unique(df, :userid))[1]
    vrand = rand(0:100, 1, length(users))

    for i=1:length(users)
        if p.per_infected  <= vrand[i]
            vstatus[i] = 0
        end
    end

    push!(
        per_infected_data,
        p.per_infected => vec(vstatus)
    )
end


#########################
# Simulation
########################
simulation_data = Dict{String, Array{Pair{String, NamedTuple}, 1}}()

infectedXloc_data = Dict{String, Any}()

for testtype in keys(test_data)
    for test in get(test_data, testtype, nothing)
        to_print = string(
            "Experiment code = $(test[:exp_id]) | Configuration label = $(test[:label]) | Perc infected = $(test[:per_infected]) | ",
            "Δ = $(test[:Δ]) | δ = $(test[:δ]) | βd = $(test[:βd]) | βᵢ = $(test[:βᵢ]) | βₑ = $(test[:βₑ]) | γₑ = $(test[:γₑ]) | γₐ = $(test[:γₐ])"
        )
        println(to_print)

        runningparams = get(intervals_data, "$(test[:Δ])$(test[:δ])", Dict{Symbol, Any}())

        SIS_per_infected_sim, infectedXloc =
            simulate(
                TVHEpidemicDynamics.SIS_ble(),
                get!(runningparams, :df, nothing),
                get!(runningparams, :intervals, nothing),
                get!(runningparams, :user2vertex, nothing),
                get!(runningparams, :loc2he, nothing),
                convert(Dates.Millisecond, Dates.Minute(test[:δ]));
                Δ = test[:Δ],
                vstatus = per_infected_data[test[:per_infected]],
                per_infected = test[:per_infected],
                c = test[:c],
                βd = test[:βd],
                βᵢ = test[:βᵢ],
                βₑ = test[:βₑ],
                γₑ = test[:γₑ],
                γₐ = test[:γₐ],
                niter = 10,
                output_path = "$(output_path)/csv/$(test[:exp_id])_$(test[:data])_$(Dates.format(now(), "Y-mm-ddTHH-MM-SS")).csv"
            )

        # get the average over all iterations
        infected_distribution = mean(collect(values(SIS_per_infected_sim)))

        push!(
            get!(simulation_data, testtype, Array{Dict{String, NamedTuple}, 1}()),
            test[:label] => (
                infected_distribution = infected_distribution,
                Δ = test[:Δ],
                δ = test[:δ],
                ncheckins = intervals_data["$(test[:Δ])$(test[:δ])"][:checkins])
            )

        push!(infectedXloc_data, testtype => infectedXloc)
    end
end


#########################
# Plotting infected ditribution
########################
linestyles = ["solid", "dashed", "dashdot", "dotted"]
markers = ["", "", "", "", "x", "+"]

for test_type in keys(simulation_data)
    linestyle = 1
    marker = 1
    labels = Array{String, 1}()
    mytitle = "$(test_type)_$(Dates.format(now(), "Y-mm-ddTHH-MM-SS")).png"

    clf()
    fig = plt.figure(figsize=(7,4))
    ax = fig.add_subplot(111)

    ncheckins = Array{Int, 1}()

    for exp in get!(simulation_data, test_type, Array{Float64, 1}())
        ylim(bottom=0.0)#, top=0.6)
        plot(exp.second.infected_distribution, linestyle=linestyles[linestyle], marker=markers[marker], markevery=10, markersize=6.5)

        xlabel("Time intervals", fontweight="semibold", labelpad=10, fontsize="x-large")
        ylabel("Δ = $(exp.second.Δ) hours \n Infected nodes in %", fontweight="semibold", fontsize="x-large", labelpad=10)
        title("δ = $(exp.second.δ) minutes", pad=10, fontweight="semibold", fontsize="x-large")

        tick_params(labelsize="large")

        push!(labels, exp.first)

        linestyle = (linestyle + 1) % (length(linestyles)+1)
        marker = (marker + 1) % (length(markers)+1)

        if linestyle == 0
            linestyle = 1
        end
        if marker == 0
            marker = 1
        end

        ncheckins = exp.second.ncheckins
    end

    legend(labels, fontsize="large", ncol=2)

    # ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis

    # color = "tab:blue"
    # ax2.set_ylabel("Number of checkins", color=color)  # we already handled the x-label with ax1
    # ax2.plot(ncheckins, color=color)
    # ax2.tick_params(axis="y", labelcolor=color)
    
    plt.tight_layout(.5)
    savefig("$(output_path)/plot/$(mytitle)")
end

gcf()



#########################
# Plotting new infected
# per location per day
########################
using Serialization
serialize("ble_data_heatmap.data", data_to_plot)

data_to_plot = infectedXloc_data["exp415"]
sim_data = intervals_data["415"]

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
savefig("heatmap.png")

