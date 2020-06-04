"""
Returns number of infected indiduals in a given location (hyperedge)
"""
function infected(h, he, vstatus, istatus)
    vertices = getvertices(h, he)
    sum = 0
    for v in vertices
        if istatus[v.first] == 0
            sum += vstatus[v.first]
        end
    end
    sum
end


"""
Regularizes the number of infected individuals
"""
function f(ninfected, c)
    ninfected > c ? c : ninfected
end


"""
Simulation process on TVH for the SIS model, considering direct and indirect contacts
between individuals and locations.
"""
function TVHSIS(
    df,
    intervals,
    node_index_map,
    he_index_map,
    Δ,
    δ;
    path = "",
    mytitle = "",
    vstatus = nothing,
    per_infected = 20,
    c = 12,
    τd = 0.2,
    τₕ = 0.2,
    τᵥ = 0.20,
    γₕ = 0.05,
    γᵥ = 0.05,
    αᵥ = 0.05,
    βᵥ = 0.0,
    αₑ = 0.05,
    βₑ = 0.0,
    ntrials = 1,
    immunizenode = false,
    immunizeedge = false,
    node_imm_strategy,
    he_imm_strategy,
    )

    per_infected_sim = Dict{Int,Array{Float64,1}}()
    per_infected_sim_trials = Dict{Int,Array{Float64,1}}()

    fname = isnothing(c) ? "AAMAS/results/1year-Δ:$(Δ)-perinfected:$(per_infected)-c:median-τₕ:$(τₕ)-τᵥ:$(τᵥ)-γₕ:$(γₕ)-γᵥ:$(γᵥ)-αᵥ:$(αᵥ)-βᵥ:$(βᵥ).csv" : "$(path)$(Base.Filesystem.path_separator)$(mytitle).csv"

    open(fname, "a") do fhandle
        write(
        fhandle,
        "time,delta,newusers,movedusers,c,td,th,tv,gh,gv,aq,bq,aeq,beq,density
        ,avg_degree,perc_tot_infected,tot_he_infected,tot_he_infected_1\n")

        for trial = 1:ntrials
            println(trial)

            # random status initialization
            h = nothing
            added, moved = 0, 0

            numnodes = length(keys(node_index_map))
            numhe = length(keys(he_index_map))
            usersepoc = rand(0:0, 1, numnodes)

            if isnothing(vstatus)
                vstatus = rand(1:1, 1, numnodes)
                vrand = rand(0:100, 1, numnodes)
                for i = 1:numnodes
                    if per_infected  <= vrand[i]
                        vstatus[i] = 0
                    end
                end
            end
            hstatus = rand(0:0, 1, numhe)

            # immunization
            istatus = rand(0:0, 1, numnodes)
            irand = rand(0:100, 1, numnodes)

            ihestatus = rand(0:0, 1, numhe)
            iherand = rand(0:100, 1, numhe)

            vnextstatus = zeros(Int64, 1, numnodes)
            nextistatus = zeros(Int64, 1, numnodes)

            henextstatus = zeros(Int64, 1, numhe)
            nextihestatus = zeros(Int64, 1, numhe)

            push!(
                get!(per_infected_sim, 0, Array{Float64,1}()),
                sum(vstatus) / length(vstatus))

            push!(
                get!(per_infected_sim_trials, 0, Array{Float64,1}()),
                sum(vstatus) / length(vstatus))

            # simulate
            for t = 1:length(intervals)
                h, a, m = generatehg(h, df, get(intervals, t, 0).first, get(intervals, t, 0).second, node_index_map, he_index_map, t, true, usersepoc)

                if t == 200
                    if immunizenode
                        for i = 1:numnodes
                            if αᵥ > irand[i] / 100
                                nextistatus[i] = 1
                            end
                        end
                    elseif immunizeedge
                        for i = 1:numhe
                            if αₑ > iherand[i] / 100
                                nextihestatus[i] = 1
                            end
                        end
                    end
                end

                if isnothing(h)
                    continue
                end

                dist = Array{Int,1}()
                if isnothing(c)
                    for he = 1:size(h)[2]
                        push!(dist, length(getvertices(h, he)))
                    end
                    c = median(dist)
                    println(t, " -- ", c)
                end

                # number of possibile interactions in the hypergraph
                density = 0
                for he = 1:size(h)[2]
                    density += length(getvertices(h, he))
                end

                # hg avg node degree
                avgdegree = 0
                for v = 1:nhv(h)
                    avgdegree += length(gethyperedges(h, v))
                end
                avgdegree /= length(node_index_map)

                hstatus_greater_1 = 0
                for he = 1:size(h)[2]
                    if hstatus[he] == 1 && length(getvertices(h, he)) > 1
                        hstatus_greater_1 += 1
                    end
                end

                for iter = 1:1
                    # direct influence
                    avgdirectdegree = 0
                    for v = 1:size(h)[1]
                        node_imm_strategy(v) && continue
                        if usersepoc[v] == 1

                            i = 0
                            for he in gethyperedges(h, v)
                                for u in getvertices(h, he.first)
                                    if v != u.first
                                        if abs(h[v, he.first] - h[u.first, he.first]) <= δ.value
                                            if vstatus[v] == 0
                                                i += vstatus[u.first]
                                            end
                                            avgdirectdegree += 1
                                        end
                                    end
                                end
                            end
                            if vstatus[v] == 0 && rand(1)[1] < 1 - ℯ^-(τd * i)
                                vnextstatus[v] = 1
                            end

                        end
                    end

                    avgdirectdegree \= sum(usersepoc)

                    #
                    # undirect influence
                    #
                    for he = 1:size(h)[2]
                        he_imm_strategy(he) && continue

                        if ihestatus[he] == 1
                            continue
                        end
                        if length(getvertices(h, he)) > 1 && hstatus[he] == 0
                            i = infected(h, he, vstatus, istatus)
                            if rand(1)[1] <  1 - ℯ^- (τₕ * f(i, c))
                                henextstatus[he] = 1
                            end
                        elseif rand(1)[1] <  1 - ℯ^- γₕ
                            henextstatus[he] = 0
                        end
                    end

                    for v = 1:size(h)[1]
                        if istatus[v] == 1 # if individual is immunized
                            continue
                        end
                        if usersepoc[v] == 1
                            if vstatus[v] == 0
                                i = 0
                                for he in gethyperedges(h, v)
                                    if length(getvertices(h, he.first)) > 1
                                        if ihestatus[he.first] == 0
                                            i += hstatus[he.first]
                                        end
                                    end
                                end
                                if rand(1)[1] < 1 - ℯ^-(τᵥ * f(i, c))
                                    vnextstatus[v] = 1
                                end
                            elseif rand(1)[1] < 1 - ℯ^- γᵥ # healthy and out of quarentine
                                vnextstatus[v] = 0
                            end
                        end
                    end
                    vstatus = copy(vnextstatus)
                    istatus = copy(nextistatus)

                    hstatus = copy(henextstatus)
                    ihestatus = copy(nextihestatus)

                    push!(
                        get!(per_infected_sim, iter + (t - 1) * 1, Array{Float64,1}()),
                        sum(vstatus) / (length(vstatus) - sum(istatus))
                    )

                    write(fhandle, "$(t),$(Δ),$(a),$(m),$(c),$(τₕ),$(τᵥ),$(γₕ),$(γᵥ),$(αᵥ),$(βᵥ),$(αₑ),$(βₑ),$(density),$(avgdegree),$(sum(vstatus) / length(vstatus)),$(sum(hstatus)),$(hstatus_greater_1)\n")
                end
            end
        end
    end
    per_infected_sim# , hgs
end
