using DataFrames
using CSV

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

d_cat = Dict([(i, count(x->x==i, cat)) for i in u_cat])

s_cat = sort(collect(d_cat), rev=true, by=x->x[2])

CSV.write("categories.csv", s_cat)
