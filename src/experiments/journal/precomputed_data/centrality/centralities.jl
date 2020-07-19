using Pkg
Pkg.activate(".")
using TVHEpidemicDynamics
using LightGraphs
using SimpleHypergraphs
using Dates

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

usersepoc = zeros(Int, length(user2vertex))
h = nothing

h, added, moved = generatehg!(
                    h,
                    df,
                    get(intervals, 30, 0).first,
                    get(intervals, 30, 0).second,
                    user2vertex,
                    loc2he,
                    usersepoc,
                    30
                )

m = SimpleHypergraphs.adjacency_matrix(h; s=2)

g = LightGraphs.SimpleGraph(m)
bc = LightGraphs.betweenness_centrality(g)

to_return = sort!(collect(1:nhv(h)), by = x -> bc[x], rev = true)

using Serialization

serialize(to_return, "bc2.data")
