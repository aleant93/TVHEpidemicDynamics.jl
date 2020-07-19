using TVHEpidemicDynamics
using CSV
using DataFrames
using Dates
using SimpleHypergraphs

#DATASET - FOURSQUARE
dataset = "data/dataset_TSMC2014_TKY.txt"#"data/dataset_TSMC2014_TKY.txt"
header = [:userid, :venueid, :catid, :catname, :lat, :lng, :timezoneoffset, :UTCtime]
dateformat = "e u d H:M:S +0000 Y"

df = CSV.read(
    dataset;
    copycols = true,
    header = header,
    dateformat = dateformat
)

df = dropmissing(df, :UTCtime)

cat = df[!, :catname]
u_cat = unique(cat)

# FILTER
mindate = Dates.DateTime("2012-05-07T00:00:00")
maxdate = Dates.DateTime("2012-06-07T00:00:00")

df = filter(
    r-> (isnothing(mindate) || r[:UTCtime] >= mindate) &&
        (isnothing(maxdate) || r[:UTCtime] <= maxdate),
        df
)

# TAGS
tags = CSV.read("data/tags/tags.csv"; datarow=1)

social = Set(CSV.read("data/tags/social.csv"; datarow=1)[:, 1])
transportation = Set(CSV.read("data/tags/transportation.csv"; datarow=1)[:, 1])
work = Set(CSV.read("data/tags/work.csv"; datarow=1)[:, 1])


# CATEGORY, MACROCATEGORY MAPPING
# social = 1
# transportation = 2
# work = 3
cat2main = Dict{String, Int}()

for c in u_cat
    if c in social
        push!(cat2main, c => 1)
    elseif c in transportation
        push!(cat2main, c => 2)
    elseif c in work
        push!(cat2main, c => 3)
    else
        push!(cat2main, c => 0)
    end
end

l = filter(x -> x[2] == 0, cat2main)


#MAPPING
# vertices
users = unique(df, :userid)
numusers = nrow(users)

# hyperedges - upper buond
places = unique(df, :venueid)
numplaces = nrow(places)

# mapping user -> node
user2vertex = Dict{String, Int}()
index = 1

# add vertices
for user in eachrow(users)
    global index
    uid = string(user[:userid])
    if !haskey(user2vertex, uid)
        push!(user2vertex, uid => index)
        index += 1
    end
end

#BUILD HG
h = Hypergraph{Int, String, String}(length(keys(user2vertex)), 3)

for checkin in eachrow(df)

    get(cat2main, checkin.catname, nothing) == 0 && continue

    setindex!(
        h,
        1, # checkin to store
        get(user2vertex, string(checkin.userid), -1), # node id
        get(cat2main, checkin.catname, nothing) # hyperedge id
    )
end

degrees = map(v -> length(gethyperedges(h, v)), 1:nhv(h))

degrees[degrees .== 3]
degrees[degrees .== 2]
degrees[degrees .== 1]
degrees[degrees .== 0]
