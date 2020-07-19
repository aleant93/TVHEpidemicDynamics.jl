using TVHEpidemicDynamics
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
fparams = "src/experiments/AAMAS20/configs/53/aamas53.csv"
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

# Weeplace dataset
dataset = "weeplaces/weeplace_checkins_cp.csv"
header = [:userid, :venueid, :UTCtime, :lat, :lng, :city, :category]
dateformat = "yyyy-mm-ddTHH:MM:SS" #y-m-dTH:M:SZ

pth = joinpath(dirname(pathof(TVHEpidemicDynamics)),"..","data", dataset)

# The simulation will consider only the data
# within this time intervals
firstcheckindate = Dates.DateTime("2010-07-29T00:00:00")
lastcheckindate = Dates.DateTime("2010-08-28T00:00:00")

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
            pth,
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

for testtype in keys(test_data)
    for test in get(test_data, testtype, nothing)
        to_print = string(
            "Experiment code = $(test[:exp_id]) | Configuration label = $(test[:label]) | Perc infected = $(test[:per_infected]) | ",
            "Δ = $(test[:Δ]) | δ = $(test[:δ]) | βd = $(test[:βd]) | βᵢ = $(test[:βᵢ]) | βₑ = $(test[:βₑ]) | γₑ = $(test[:γₑ]) | γₐ = $(test[:γₐ])"
        )
        println(to_print)

        runningparams = get(intervals_data, "$(test[:Δ])$(test[:δ])", Dict{Symbol, Any}())

        SIS_per_infected_sim =
            simulate(
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
                niter = 1,
                output_path = "$(output_path)/csv/$(test[:exp_id])_$(test[:data])_$(Dates.format(now(), "Y-mm-ddTHH-MM-SS")).csv"
            )

        # get the average over all iterations
        infected_distribution = mean(collect(values(SIS_per_infected_sim)))

        push!(
            get!(simulation_data, testtype, Array{Dict{String, NamedTuple}, 1}()),
            test[:label] => (infected_distribution = infected_distribution, Δ = test[:Δ], δ = test[:δ])
        )
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
    figure(figsize=(7,4))

    for exp in get!(simulation_data, test_type, Array{Float64, 1}())
        ylim(bottom=0.0, top=0.6)
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
    end
    legend(labels, fontsize="large", ncol=2)
    plt.tight_layout(.5)
    savefig("$(output_path)/plot/$(mytitle)")
end

gcf()
