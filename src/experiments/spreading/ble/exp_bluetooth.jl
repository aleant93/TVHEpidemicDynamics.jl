using Pkg
Pkg.activate(".")
using TVHEpidemicDynamics
using DataFrames
using Dates
using CSV
using JSON3
using JSONTables
using PyPlot
using Statistics

"""
    Experiments on an epidemic spreading via TVH,
    also evaluating the number of checkins per time interval
"""

############################
# Loading simulation params
############################
project_path = dirname(pathof(TVHEpidemicDynamics))

# BLE BLUETOOTH
output_path = joinpath(project_path, "experiments", "spreading", "ble", "results")

fparams =
    joinpath(project_path, "experiments", "spreading", "ble", "configs", "blebluetooth.json")

fdata_params =
    joinpath(project_path, "experiments", "spreading", "ble", "configs", "blebluetooth_dataset.json")


jtable = jsontable(read(open(fparams, "r")))
paramsdf = DataFrame(jtable)

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
data_params = JSON3.read(read(open(fdata_params, "r")))
header = [Symbol(col) for col in data_params.header]

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
            data_params.dataset,
            header,
            Symbol(data_params.userid),
            Symbol(data_params.venueid),
            Symbol(data_params.UTCtime),
            data_params.dateformat;
            Δ = convert(Dates.Millisecond, Dates.Hour(i.Δ)),
            δ = convert(Dates.Millisecond, Dates.Minute(i.δ)),
            maxdate = Dates.DateTime(data_params.end_date),
            mindate = Dates.DateTime(data_params.start_date)
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
            r-> (isnothing(mindate) || r[:timestamp] >= mindate) &&
                (isnothing(maxdate) || r[:timestamp] <= maxdate),
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
per_infected = unique(paramsdf, [:infected_percentage])[!, [:infected_percentage]]
per_infected_data = Dict{Float64, Array{Int, 1}}()

users = keys(intervals_data[collect(keys(intervals_data))[1]][:user2vertex])

for p in eachrow(per_infected)
    vstatus = fill(1, length(users))
    vrand = rand(Float64, length(users))

    for i=1:length(users)
        if p.infected_percentage  <= vrand[i]
            vstatus[i] = 0
        end
    end

    push!(
        per_infected_data,
        p.infected_percentage => vec(vstatus)
    )
end


#########################
# Simulation
########################
simulation_data = Dict{String, Array{Pair{String, NamedTuple}, 1}}()

infectedXloc_data = Dict{String, Any}()

for testtype in keys(test_data)
    for test in get(test_data, testtype, nothing)

        println("----------------EXP CONFIG-------------------------")
        for property in propertynames(test)
            print("$(property) = $(test[property])  |   ")
        end
        println("\n---------------------------------------------------")

        runningparams = get(intervals_data, "$(test[:Δ])$(test[:δ])", Dict{Symbol, Any}())

        res_path =
            joinpath(output_path, "csv", "$(test[:exp_id])_$(test[:exp])_$(Dates.format(now(), "Y-mm-ddTHH-MM-SS")).csv")

        SIS_per_infected_sim, infectedXloc =
            simulate(
                SIS_infected_loc(),
                get!(runningparams, :df, nothing),
                get!(runningparams, :intervals, nothing),
                get!(runningparams, :user2vertex, nothing),
                get!(runningparams, :loc2he, nothing),
                convert(Dates.Millisecond, Dates.Minute(test[:δ]));
                Δ = test[:Δ],
                vstatus = per_infected_data[test[:infected_percentage]],
                per_infected = test[:infected_percentage],
                c = test[:c],
                βd = test[:βd],
                βᵢ = test[:βᵢ],
                βₑ = test[:βₑ],
                γₑ = test[:γₑ],
                γₐ = test[:γₐ],
                niter = 10,
                output_path = res_path
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
        ylim(bottom=0.0, top=2)
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

    ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis

    color = "tab:blue"
    ax2.set_ylabel("Number of checkins", color=color)  # we already handled the x-label with ax1
    ax2.plot(ncheckins, color=color)
    ax2.tick_params(axis="y", labelcolor=color)

    plt.tight_layout(.5)
    savefig("$(output_path)/plot/$(mytitle)")
end

gcf()


# HEATMAP DATA

# data_to_plot = infectedXloc_data["exp415-low"]
# sim_data = intervals_data["415"]
#
# serialize("ble_data_heatmap_infectedPerLocation.data", infectedXloc_data)
# serialize("sim_data_infectedPerLocation.data", intervals_data)
