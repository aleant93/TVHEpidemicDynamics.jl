using Pkg
Pkg.activate(".")
using TVHEpidemicDynamics
using DataFrames
using Dates
using CSV
using JSON3
using PyPlot
using Statistics
using StatsBase

"""
    This script evaluates:
    - which rooms to close
    - which people can enter
    implementing a policy of social distancing
"""
project_path = dirname(pathof(TVHEpidemicDynamics))

output_path = joinpath(project_path, "experiments", "immunization", "ble", "results")

# Simulation parameters
fparams =
    joinpath(project_path, "experiments", "spreading", "ble", "configs", "blebluetooth.csv")

fdata_params =
    joinpath(project_path, "experiments", "spreading", "ble", "configs", "blebluetooth_dataset.json")


data_params = JSON3.read(read(open(fdata_params, "r")))

# The choice of the interval within which
# either an indirect (Δ) or direct (δ) contact
# may occur influences the data the
# simulation is run on.
# For this reason, it is necessary to store
# diffent information according to the
# values of both Δ and δ.
intervals_data = Dict{String, Dict{Symbol, Any}}()

header = [Symbol(col) for col in data_params.header]

data_params.start_date

# evaluating new dataset
# removing people
df, intervals, user2vertex, loc2he =
        generate_model_data(
            data_params.dataset,
            header,
            Symbol(data_params.userid),
            Symbol(data_params.venueid),
            Symbol(data_params.UTCtime),
            data_params.dateformat;
            Δ = convert(Dates.Millisecond, Dates.Hour(24)),
            δ = convert(Dates.Millisecond, Dates.Minute(15)),
            maxdate = Dates.DateTime(data_params.end_date),
            mindate = Dates.DateTime(data_params.start_date)
        )

# number of peole checking-in within one day
n_people = Array{Int, 1}()

for i in 1:length(intervals)

    _df = filter(
        r-> (r[:timestamp] >= intervals[i].first) &&
            (r[:timestamp] <= intervals[i].second),
            df
    )

    users = unique(df, Symbol(data_params.userid))

    push!(n_people, nrow(users))
end

most_crowded_day = findmax(n_people)
intervals[most_crowded_day[2]]

# we will only allow 31 people
# to access the building
df_cp = copy(df)
checkins_to_remove = Set{Int}()

for i in 1:length(intervals)
    added = 0

    _df = filter(
        r-> (r[:timestamp] >= intervals[i].first) &&
            (r[:timestamp] <= intervals[i].second),
            df
    )

    users = unique(df, Symbol(data_params.userid))
    chosen = sample(users[!, :userid], 31)

    for checkin in eachrow(_df)
        if !(checkin.userid in chosen)
            push!(checkins_to_remove, checkin.entry_id)
        end

    end
end

to_store = filter(
        row -> !(row.entry_id in checkins_to_remove),
        df
    )

CSV.write("ble.csv", to_store, dateformat="Y-m-dTH:M:S")

CSV.read("ble.csv", dateformat="Y-m-dTH:M:S")

#
# CLOSING ROOMS
#
order = ["00" "01" "02" "30" "31" "25" "29" "08" "13" "28" "24" "23" "27" "26" "22" "14" "15" "16" "17" "18" "19" "20" "21" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12"]
_ids = ["rpi-$(id)" for id in order]

#rooms_to_close = _ids[1:2:length(_ids)]
rooms_to_close = ["30" "01" "29" "13" "24" "27" "22" "20" "18" "15" "5" "7" "11"]
_ids_rooms_to_close = ["rpi-$(id)" for id in rooms_to_close]

to_store_rooms = filter(
    row -> !(row.venueid in _ids_rooms_to_close),
    df
)

CSV.write("ble_closed_rooms.csv", to_store_rooms, dateformat="Y-m-dTH:M:S")


for row in eachrow(df)
    if row.venueid in rooms_to_close
end



###
# interval => people to verify whether they are infected or not
people_to_check = Dict{DateTime, Array{String, 1}}()

for i in 1:length(intervals)

    _df = filter(
        r-> (r[Symbol(data_params.UTCtime)] >= intervals[i].first) &&
            (r[Symbol(data_params.UTCtime)] <= intervals[i].second),
            df
    )

    users = unique(df, Symbol(data_params.userid))[!, Symbol(data_params.userid)]
    #chosen = sample(users[!, :userid], 10)
    users_string = [string(u) for u in users]

    push!(people_to_check, intervals[i].first => users_string)
end

using Serialization
serialize("people_to_check.data", people_to_check)

_df = filter(
    r-> (r[Symbol(data_params.UTCtime)] >= intervals[1].first) &&
        (r[Symbol(data_params.UTCtime)] <= intervals[1].second),
        df
)

user_groups = groupby(_df, Symbol(data_params.userid))

data = sort!(_df, [order(:userid), :UTCtime])
data[data[:userid] .== 0, :]
