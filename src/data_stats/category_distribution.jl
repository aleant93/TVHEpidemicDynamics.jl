using DataFrames
using CSV

############################
# Loading data set params
############################
project_path = dirname(pathof(TVHEpidemicDynamics))

# FOURSQUARE
output_path = joinpath(project_path, "data_stats", "AAMAS20", "plots")
fparams = joinpath(project_path, "data_stats", "AAMAS20", "configs", "foursquare.json")

# WEEPLACES
output_path = joinpath(project_path, "data_stats", "weeplaces", "plots")
fparams = joinpath(project_path, "data_stats", "weeplaces", "configs", "weeplaces.json")



############################
# Reading data
############################
data_params = JSON3.read(read(open(fparams, "r")))
header = [Symbol(col) for col in data_params.header]

df = CSV.read(
    data_params.dataset;
    copycols = true,
    header = header,
    dateformat = data_params.dateformat
    )

df = dropmissing(df)

cat = df[!, :catname]
u_cat = unique(cat)

d_cat = Dict([(i, count(x->x==i, cat)) for i in u_cat])

s_cat = sort(collect(d_cat), rev=true, by=x->x[2])

CSV.write("categories.csv", s_cat)
