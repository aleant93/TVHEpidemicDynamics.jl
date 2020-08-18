using Pkg
Pkg.activate(".")
using CSV
using TVHEpidemicDynamics
using SimpleHypergraphs
using Dates
using JSON3
using StringEncodings
using Serialization

############################
# Loading data set params
# foursquare
############################
project_path = dirname(pathof(TVHEpidemicDynamics))

fdata_params =
    joinpath(project_path, "experiments", "spreading", "AAMAS20", "configs", "foursquare.json")

output_path =
    joinpath(project_path, "experiments", "immunization", "foursquare", "precomputed_data", "lockdown")

data_params = JSON3.read(read(open(fdata_params, "r")))
header = [Symbol(col) for col in data_params.header]


Δ = 4
δ = 1

df, intervals, user2vertex, loc2he =
    generate_model_data(
        data_params.dataset,
        header,
        Symbol(data_params.userid),
        Symbol(data_params.venueid),
        Symbol(data_params.UTCtime),
        data_params.dateformat;
        Δ = convert(Dates.Millisecond, Dates.Hour(Δ)),
        δ = convert(Dates.Millisecond, Dates.Minute(δ)),
        maxdate = Dates.DateTime(data_params.end_date),
        mindate = Dates.DateTime(data_params.start_date)
    )

size(df)

df = unique(df, :venueid)
df_cat = df[!, [:venueid, :catname]]

# unique_cat = unique(df, :catname)[!, [:catname]]
# CSV.write("data/tags/all_categories.csv", unique_cat)

# TAGS
tags = Set(CSV.read("data/tags/all_categories.csv"; datarow=1)[:, 1])

social = Set(CSV.read("data/tags/social.csv"; datarow=1)[:, 1])
transportation = Set(CSV.read("data/tags/transportation.csv"; datarow=1)[:, 1])
work = Set(CSV.read("data/tags/work.csv"; datarow=1)[:, 1])


social_venue = Array{Int, 1}()
transp_venue = Array{Int, 1}()
work_venue = Array{Int, 1}()
all_tags = Array{Int, 1}()

# SAVING LOC ID TO IMMUNIZE
for row in eachrow(df_cat)
    if row.catname in social || startswith(row.catname, "Caf")
        push!(social_venue, get(loc2he, row.venueid, nothing))
    elseif row.catname in transportation
        push!(transp_venue, get(loc2he, row.venueid, nothing))
    elseif row.catname in work
        push!(work_venue, get(loc2he, row.venueid, nothing))
    else
        println("$(row.venueid) -- $(row.catname)")
    end

    push!(all_tags, get(loc2he, row.venueid, nothing))
end

all_venues = vcat(social_venue, work_venue, transp_venue)

serialize(joinpath(output_path, "social_venues.data"), social_venue)
serialize(joinpath(output_path, "transportation_venues.data"), transp_venue)
serialize(joinpath(output_path, "work_venues.data"), work_venue)
serialize(joinpath(output_path, "union.data"), all_venues)
serialize(joinpath(output_path, "all_venues.data"), all_tags)

deserialize(joinpath(output_path, "social_venues.data"))
deserialize(joinpath(output_path, "transportation_venues.data"))
deserialize(joinpath(output_path, "work_venues.data"))
