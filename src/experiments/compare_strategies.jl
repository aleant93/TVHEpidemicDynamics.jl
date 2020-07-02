using TVHEpidemicDynamics
using Dates
using CSV
using PyPlot
using Statistics

#########################
# Loading simulation params
########################

# Simulation parameters
fparams = "src/experiments/configs/epidemic_control/cmp_strat.csv"

paramsdf = CSV.read(
                fparams;
                copycols=true,
                header=[:exp_id, :data, :label, :per_infected, :c,  :Δ, :δ ,:βd, :βᵢ, :βₑ, :γₑ, :γₐ, :αₑ],
                datarow=2
            )

# just a trick to group together
# all experiments to show in the same plot
test_data = Dict{String,Array{Any,1}}()
for params in eachrow(paramsdf)
    push!(
            get!(test_data, params[:exp_id], Array{Any,1}()),
            params
        )
end


#########################
# Generating model data
########################

# Foursqaure dataset 
dataset = "/home/antonio/Scrivania/Tesi/data/dataset_TSMC2014_TKY.txt"
header = [:userid, :venueid, :catid, :catname, :lat, :lng, :timezoneoffset, :UTCtime]
dateformat = "e u d H:M:S +0000 Y"

# The simulation will consider only the data
# within this time intervals
firstcheckindate = Dates.DateTime("2012-05-07T00:00:00")
lastcheckindate = Dates.DateTime("2012-06-07T00:00:00")

Δₕ = 12
Δₘ = convert(Dates.Millisecond, Dates.Hour(Δₕ))

δₛ = 60
δₘ = convert(Dates.Millisecond, Dates.Minute(δₛ))

df, intervals, node_index_map, he_index_map =
        generate_model_data(
            dataset, header, dateformat;
            Δ=Δₘ, δ=δₘ,
            mindate=firstcheckindate, maxdate=lastcheckindate)

#########################
# Simulation
########################
simulation_data = Dict{String,Array{Pair{String,Array{Float64,1}},1}}()

for testtype in keys(test_data)
    for test in get(test_data, testtype, nothing)
        to_print = string(
                "Experiment code = $(test[:exp_id]) | Configuration label = $(test[:label]) | Perc infected = $(test[:per_infected]) | ",
                "Δ = $(test[:Δ]) | δ = $(test[:δ]) | βd = $(test[:βd]) | βᵢ = $(test[:βᵢ]) | βₑ = $(test[:βₑ]) | γₑ = $(test[:γₑ]) | γₐ = $(test[:γₐ])",
                "| αₑ = $(test[:αₑ])"
            )
        println(to_print)

        # Execute the simulation first without immunization
        SIS_per_infected_sim =
                    TVHSIS(
                        df,
                        intervals,
                        node_index_map,
                        he_index_map,
                        convert(Dates.Millisecond, Dates.Minute(test[:δ]));
                        Δ=test[:Δ],
                        c=5, 
                        βd=test[:βd],
                        βᵢ=test[:βᵢ], 
                        βₑ=test[:βₑ], 
                        γₑ=test[:γₑ], 
                        γₐ=test[:γₐ],
                        αₑ=test[:αₑ],
                        niter=1,
                        output_path="src/experiments/results/epidemic_control/nodes_vs_hes/partial/csv/$(test[:exp_id])_$(test[:data])_$(Dates.format(now(), "Y-mm-ddTHH-MM-SS")).csv"
                        )

        # get the average over all iterations 
        infected_distribution = mean(collect(values(SIS_per_infected_sim)))

        push!(
            get!(simulation_data, testtype, Array{Dict{String,Array{Float64,1}},1}()),
            "no control" => infected_distribution
        ) 

        # Then apply each immunization strategy to nodes and hes
        for target in ["nodes", "hes"]
            for strategy in [uniform, contacts_based, acquaintance]
                nodes_imm_strategy = target == "nodes" ? strategy : nothing
                hes_imm_strategy = target == "hes" ? strategy : nothing

                SIS_per_infected_sim =
                    TVHSIS(
                        df,
                        intervals,
                        node_index_map,
                        he_index_map,
                        convert(Dates.Millisecond, Dates.Minute(test[:δ]));
                        Δ=test[:Δ],
                        c=5,
                        βd=test[:βd],
                        βᵢ=test[:βᵢ],
                        βₑ=test[:βₑ],
                        γₑ=test[:γₑ],
                        γₐ=test[:γₐ],
                        αₑ=test[:αₑ],
                        niter=1,
                        output_path="src/experiments/results/epidemic_control/nodes_vs_hes/partial/csv/$(test[:exp_id])_$(test[:data])_$(Dates.format(now(), "Y-mm-ddTHH-MM-SS")).csv",
                        nodes_imm_strategy=nodes_imm_strategy,
                        hes_imm_strategy=hes_imm_strategy
                        )

            # get the average over all iterations 
                infected_distribution = mean(collect(values(SIS_per_infected_sim)))

                push!(
                get!(simulation_data, testtype, Array{Dict{String,Array{Float64,1}},1}()),
                target * " " * string(strategy) => infected_distribution
            ) 
            end
        end
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
    labels = Array{String,1}()
    mytitle = "$(test_type)_$(Dates.format(now(), "Y-mm-ddTHH-MM-SS")).png"

    clf()
    figure(figsize=(7, 4))
    
    for exp in get!(simulation_data, test_type, Array{Float64,1}())        
        ylim(bottom=0.0, top=0.6)
        plot(exp.second, linestyle=linestyles[linestyle], marker=markers[marker], markevery=10, markersize=6.5)

        ylabel("Infected nodes in %", fontweight="semibold", fontsize="x-large", labelpad=10)
        xlabel("Time intervals", fontweight="semibold", labelpad=10, fontsize="x-large")

        tick_params(labelsize="large")

        push!(labels, exp.first)

        linestyle = (linestyle + 1) % (length(linestyles) + 1)
        marker = (marker + 1) % (length(markers) + 1)

        if linestyle == 0
            linestyle = 1
        end
        if marker == 0
            marker = 1
        end
    end
    legend(labels, fontsize="large", ncol=2)
    plt.tight_layout(.5)
    savefig("src/experiments/results/epidemic_control/nodes_vs_hes/partial/plot/$(mytitle)")
end

gcf()