using CSV
using Dates
using DataFrames
using TVHEpidemicDynamics

df = CSV.read(
    "data/blebeacon/BLEBeacon-Dataset-master/Check-In Check-Out Report.csv",
    dateformat = "Y-m-d H:M:S"
    )

dataset = "data/blebeacon/BLEBeacon-Dataset-master/Check-In Check-Out Report.csv"#"data/dataset_TSMC2014_TKY.txt"
header = [:entry_id, :beacon_id, :in_time, :out_time, :rpi_id]
dateformat = "Y-m-d H:M:S"

# The simulation will consider only the data
# within this time intervals
firstcheckindate = Dates.DateTime("2016-09-06T00:00:00")
lastcheckindate = Dates.DateTime("2016-10-19T00:00:00")

Δ = 4
δ = 1


df = CSV.read(
    dataset;
    copycols = true,
    header = header,
    dateformat = dateformat
)

df, intervals, user2vertex, loc2he =
    generate_model_data(
        dataset,
        header,
        :beacon_id,
        :rpi_id,
        :in_time,
        dateformat;
        Δ = convert(Dates.Millisecond, Dates.Hour(Δ)),
        δ = convert(Dates.Millisecond, Dates.Minute(δ)),
        maxdate = lastcheckindate,
        mindate = firstcheckindate
    )