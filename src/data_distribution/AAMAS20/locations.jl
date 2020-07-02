using TVHEpidemicDynamics
using DataFrames
using CSV
using Dates
using PyPlot
using Statistics

"""
 Studying location distribution and, consequently,
 indirect contacts distribution.
"""

dataset = "data/dataset_TSMC2014_TKY.txt"
header = [:userid, :venueid, :catid, :catname, :lat, :lng, :timezoneoffset, :UTCtime]
dateformat = "e u d H:M:S +0000 Y" #y-m-dTH:M:SZ

df = CSV.read(
    dataset;
    copycols = true,
    header = header,
    dateformat = dateformat
    )

#getting first checkin within the data
mindate = minimum(df[!, :UTCtime])

#evaluating data starting from the first monday
#minMonday = 2012-04-09T18:17:18
#minmonday = Dates.DateTime("2012-04-09T18:17:18")
minmonday = Dates.DateTime("2012-04-09T00:00:00")
mintime = minmonday


#
# DISTRIBUTION PARAMETERS
#
Δₕ = 24*7
Δₘ = convert(Dates.Millisecond, Dates.Hour(Δₕ))

# Evaluating the total number of time slots
# within the observation period (last_checkin_date - minmonday)
# using Δₘ as discretization param
nintervals = ceil(Int,(maximum(df[!, :UTCtime]) - minmonday)/Δₘ)

# interval_id -> (start_date, end_date)
intervals = Dict{Int, Pair{DateTime, DateTime}}()
evaluateintervals!(mintime, maximum(df[!, :UTCtime]), intervals, nintervals, Δₘ)


# Count checkin data according to the discrete time intervals
# evaluated at the previous step
checkinsperinterval = Dict{Int, Int}()
evaluatedensity!(intervals, checkinsperinterval, df)

# Sort time intervals according to the number of checkins
# within each interval
sortedcheckinsperinterval = sort(collect(checkinsperinterval), by=x->x[2], rev=true)

# Working with the week containing
# the majority of checkins
maxweek = get(intervals, sortedcheckinsperinterval[1].first, nothing)
#maxweek = intervals[5]


"""
 Studying location distribution
 within the most crowded week.

 Here, we will use a discretization interval of δₕ = 4 hours.
"""
δₕ = 4
δₘ = convert(Dates.Millisecond, Dates.Hour(δₕ))

# Evaluating the total number of time slots
# within the week presenting the majority of checkins
# using δₘ as discretization parameter
nintervals_δ = ceil(Int,(maxweek.second - maxweek.first)/δₘ)

# interval_id -> (start_date, end_date)
intervals_δ = Dict{Int, Pair{DateTime, DateTime}}()
evaluateintervals!(maxweek.first, maxweek.second, intervals_δ, nintervals_δ, δₘ)

# Evaluatin location distribution within each interval.
distancewithinintervals = Dict{Int, Array{Any, 1}}()
evaluate_location_distribution!(intervals_δ, distancewithinintervals, df, convert(Dates.Millisecond, Dates.Hour(1)))



"""
    Plot location distribution within 7th-13th May.
"""
clf()

fig = plt.figure(figsize=(10, 4))
ax = fig.add_subplot(111)

ylabel("Time in Seconds", fontweight="semibold")
xlabel("Time Intervals", fontweight="semibold")

for t=1:length(keys(distancewithinintervals))
    values = get(distancewithinintervals, t, Array{Any, 1}())
    ax.boxplot(values, widths=0.5, positions=[t])
end

ax.set_xlim(0, 43)

xlabels = ["7 May", "04", "08", "12", "16", "20", "8 May", "04", "08", "12", "16", "20", "9 May", "04", "08", "12", "16", "20", "10 May", "04", "08", "12", "16", "20", "11 May", "04", "08", "12", "16", "20", "12 May", "04", "08", "12", "16", "20", "13 May", "04", "08", "12", "16", "20"]
#xlabels = ["7\nMay", "", "", "12h", "", "", "8\nMay", "", "", "12h", "", "", "9\nMay", "", "", "", "", "", "10\nMay", "", "", "", "", "", "11\nMay", "", "", "", "", "", "12\nMay", "", "", "", "", "", "13\nMay", "", "", "", "", ""]

xticks(1:42, xlabels, rotation=80)

gcf()

plt.tight_layout(.5)

savefig("src/data_distribution/plots/location_distribution.png")



#
# plotting tests
#
toplot = Vector{Array{Any,1}}()
cut = 0

for t=1:length(keys(distancewithinintervals))
    values = get(distancewithinintervals, t, Array{Any, 1}())
    #println(t, "-->", values)
    #cut = 60
    # if length(values) == 0 || t==16
    #     cut = 60
    # else
    #     cut = quantile(values, 0.9)
    # end

    # try
    #     println("$(t) -- $(cut) -- $(median(values)) -- $(length(values[values.<cut])/length(values)) -- $(length(values[values.>cut])/length(values))")
    # catch ArgumentError
    #     println("$(t) -- 60 -- NaN -- $(length(values[values.<cut])/length(values)) -- $(length(values[values.>cut])/length(values))")
    # end

    #filter!(x -> x < cut, values)
    push!(toplot, values)
end


clf()
boxplot(toplot)
gcf()
savefig("src/AAMAS/plots/res/indirectcontacts.png")


#
#
#
uniqueitems = Dict{Int, Dict{Int, Int}}()

for t=1:length(keys(distancewithinintervals))
    values = get(distancewithinintervals, t, Array{Any, 1}())
    u = unique(values)
    d = Dict([(i, count(x -> x == i, values)) for i in u])

    get!(uniqueitems, t, Dict{Int, Int}())

    push!(
        uniqueitems,
        t => d
    )
end

for elem in uniqueitems
    println(elem)
end
