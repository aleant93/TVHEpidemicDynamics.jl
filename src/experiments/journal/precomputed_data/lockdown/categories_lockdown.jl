using Pkg
Pkg.activate(".")
using CSV
using TVHEpidemicDynamics
using LightGraphs
using SimpleHypergraphs
using Dates
using StringEncodings
using Serialization

# Foursqaure dataset
dataset = "data/dataset_TSMC2014_TKY.txt"#"data/dataset_TSMC2014_TKY.txt"
header = [:userid, :venueid, :catid, :catname, :lat, :lng, :timezoneoffset, :UTCtime]
dateformat = "e u d H:M:S +0000 Y"

# The simulation will consider only the data
# within this time intervals
firstcheckindate = Dates.DateTime("2012-05-07T00:00:00")
lastcheckindate = Dates.DateTime("2012-06-07T00:00:00")

Δ = 4
δ = 1

df, intervals, user2vertex, loc2he =
    generate_model_data(
        dataset,
        header,
        :userid,
        :venueid,
        :UTCtime,
        dateformat;
        Δ = convert(Dates.Millisecond, Dates.Hour(Δ)),
        δ = convert(Dates.Millisecond, Dates.Minute(δ)),
        maxdate = lastcheckindate,
        mindate = firstcheckindate
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

serialize("src/experiments/journal/precomputed_data/social_venues.data", social_venue)
serialize("src/experiments/journal/precomputed_data/transportation_venues.data", transp_venue)
serialize("src/experiments/journal/precomputed_data/work_venues.data", work_venue)
serialize("src/experiments/journal/precomputed_data/union.data", all_venues)
serialize("src/experiments/journal/precomputed_data/all_venues.data", all_tags)

deserialize("src/experiments/journal/precomputed_data/social_venues.data")
deserialize("src/experiments/journal/precomputed_data/transportation_venues.data")
deserialize("src/experiments/journal/precomputed_data/work_venues.data")
