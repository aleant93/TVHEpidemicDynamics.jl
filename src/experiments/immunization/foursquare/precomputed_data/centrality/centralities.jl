using Pkg
Pkg.activate(".")
using TVHEpidemicDynamics
using JSON3
using SimpleHypergraphs
using Dates

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
