"""
    generate_model_data()

Load a dataframe from a file `fname` with columns `columns`, containing checkin-data.
Return:
- a `dataframe` containing all the checkins that will be used for the simulation;
- a `dict` of intervals within which an indirect contact may occur:
    - at each interval `i` corresponds a pair of dates delimiting that interval
    - `i` => (date₁, date₂)
- a `dict` maping each user to the corresponding node id in the hypergraph
- a `dict` maping each location to the corresponding hyperedge id in the hypergraph

### Arguments
- `dateformat`: the format of the timestamp of the checkin;
- `Δ`: time within which an indirect contact may happen;
- `δ`: time within which a direct contact may happen;
- `limit`: number of rows to read. Each row represents a checkin;
- `mindate`: if set, a checkin is considered only if it happend after `mindate`;
- `maxdate`: if set, a checkin is considered only if it happend before `maxdate`.
"""
function generate_model_data(
    fname::AbstractString,
    columns::AbstractVector,
    dateformat::AbstractString;
    Δ::Millisecond=Dates.Millisecond(14400000),
    δ::Millisecond=Dates.Millisecond(1000),
    limit::Int=typemax(Int),
    mindate::Union{DateTime,Nothing}=nothing,
    maxdate::Union{DateTime,Nothing}=nothing
)

    df = CSV.read(
        fname;
        copycols=true,
        header=columns,
        dateformat=dateformat,
        limit=limit
    )

    # remove all rows without a timestamp
    df = dropmissing(df, :UTCtime)

    # filter out all the rows not in the 
    # timeframe to consider
    df = filter(
        r->(isnothing(mindate) || r.UTCtime >= mindate) &&
            (isnothing(maxdate) || r.UTCtime <= maxdate), 
            df
    )

    sort!(df, [:UTCtime])

    # vertices
    users = unique(df, :userid)
    numusers = nrow(users)

    # hyperedges - upper buond
    places = unique(df, :venueid)
    numplaces = nrow(places)

    # time slots
    mintime = minimum(df[!, :UTCtime])
    maxtime = maximum(df[!, :UTCtime])

    nintervals = ceil(Int, (maxtime - mintime) / Δ)
    intervals = Dict{Int,Pair{DateTime,DateTime}}()

    # build timeframes
    for i = 1:nintervals
        offset = 
            mintime + Δ + δ > maxtime ? maxtime + Dates.Millisecond(1) : mintime + Δ + δ

        get!(
            intervals, 
            i, 
            Pair{DateTime,DateTime}(
                convert(Dates.DateTime, (mintime - δ)), convert(Dates.DateTime, (offset))
            )
        )

        mintime = offset - δ
    end

    # mapping user -> node
    node_index_map = Dict{String,Int}()
    index = 1

    # add vertices
    for user in eachrow(users)
        userid = string(user[:userid])
        if !haskey(node_index_map, userid)
            push!(node_index_map, userid => index)
            index += 1
        end
    end

    # mapping place -> he
    he_index_map = Dict{String,Int}()
    index = 1

    tags = Array{String,1}()
    # add he
    for place in eachrow(places)
        placeid = string(place[:venueid])
        if !haskey(he_index_map, placeid)
            push!(he_index_map, placeid => index)
            push!(tags, place[:catname])
            index += 1
        end
    end

    df, intervals, node_index_map, he_index_map, tags
end