using Pkg
Pkg.activate(".")
using TVHEpidemicDynamics
using CSV
using DataFrames
using Dates
using JSON3
using JSONTables
using PyPlot
using Statistics
using Serialization

"""
    Experiments on the efficacy of using instant swabs
"""

############################
# Loading simulation params
############################
project_path = dirname(pathof(TVHEpidemicDynamics))

# BLE BLUETOOTH
output_path = joinpath(project_path, "experiments", "immunization", "ble", "results")

fparams =
    joinpath(project_path, "experiments", "immunization", "ble", "configs", "blebluetooth_swabs_3d.json")

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

    vstatus = fill(0, length(users))
    vstatus[1] = 1

    push!(
        per_infected_data,
        p.infected_percentage => vec(vstatus)
    )
end


#########################
# Simulation
########################
simulation_data = Dict{String, Array{Pair{String, NamedTuple}, 1}}()
people_to_check = deserialize(joinpath("data", "blebeacon", "people_to_check.data"))

swabs_data = Dict{Int, Array{Float64, 1}}()

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

        # create dict with other params

        SIS_per_infected_sim =
            simulate(
                SIS_swabs(),
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
                output_path = res_path,
                Dict{}(:n_swabs => test[:n_swabs], :people_per_day => people_to_check)...
            )

        # get the average over all iterations
        infected_distribution = mean(collect(values(SIS_per_infected_sim)))

        push!(swabs_data, test[:n_swabs] => infected_distribution)

        push!(
            get!(simulation_data, testtype, Array{Dict{String, NamedTuple}, 1}()),
            test[:label] => (infected_distribution = infected_distribution, Δ = test[:Δ], δ = test[:δ])
        )
    end
end


#########################
# Plotting a surface
# with percentage of infected users
# per timeframe
# based on the number of swabs used
########################
using Plots;
pyplot()

# x = range(-2,stop=2,length=100)
# y = range(sqrt(2),stop=2,length=100)
#
# f(x,y) = x*y-x-y+1

x = range(1, stop = length(swabs_data[0]))
y = range(0, stop = length(swabs_data) - 1)

f_infected(x, y) = swabs_data[y][x]

Plots.plot(x, y, f_infected, st=:surface, camera=(200, 30))

Plots.savefig(joinpath(output_path, "plot", "swabs.png"))
