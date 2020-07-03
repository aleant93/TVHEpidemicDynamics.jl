"""
    infected(h, he, vstatus, istatus)

Count the number of infected nodes within an hyperedge. 
"""
function infected(h, he, vstatus, istatus)
    vertices = getvertices(h, he)
    sum = 0
    for v in vertices
        # if the node is not immunized
        if istatus[v.first] == 0
            sum += vstatus[v.first]
        end
    end
    sum
end


"""
    f(num_infected, c)

Non-linear function used to bound the infection pressure 
for large values of 𝑛.
"""
function f(num_infected, c)
    num_infected > c ? c : num_infected
end


"""
TVHSIS(
    df::DataFrame,
    intervals::Dict{Int,Pair{DateTime,DateTime}},
    node_index_map::Dict{String,Int},
    he_index_map::Dict{String,Int},
    δ::Dates.Millisecond;
    Δ::Union{Int,TimePeriod,Nothing} = nothing,
    vstatus::Union{Array{Int,1},Nothing} = nothing,
    per_infected::Int = 20,
    c::Union{Int,Nothing} = 5,
    βd::Float64 = 0.2,
    βₑ::Float64 = 0.06,
    βᵢ::Float64 = 0.1,
    γₑ::Float64 = 0.06,
    γₐ::Float64 = 0.1,
    niter::Int = 1,
    output_path::Union{AbstractString,Nothing} = nothing,
    α::Float64ₑ = 0.0,
    nodes_imm_strategy=nothing,
    hes_imm_strategy=nothing,
    app_sim::Bool=false
)

Simulation process on TVH for the SIS model, considering direct and indirect contacts
between individuals and locations.
"""
function TVHSIS(
        df::DataFrame,
        intervals::Dict{Int,Pair{DateTime,DateTime}},
        node_index_map::Dict{String,Int},
        he_index_map::Dict{String,Int},
        δ::Dates.Millisecond;
        Δ::Union{Int,TimePeriod,Nothing}=nothing,
        vstatus::Union{Array{Int,1},Nothing}=nothing,
        per_infected::Int=20,
        c::Union{Int,Nothing}=5,
        βd::Float64=0.2,
        βₑ::Float64=0.06,
        βᵢ::Float64=0.1,
        γₑ::Float64=0.06,
        γₐ::Float64=0.1,
        niter::Int=1,
        output_path::Union{AbstractString,Nothing}=nothing,
        αₑ::Float64=0.0,
        nodes_imm_strategy=nothing,
        hes_imm_strategy=nothing,
        app_sim::Bool=false
)

    if isnothing(output_path)
        if !isdir("results/")
            Filesystem.mkdir("results/")
        end

        output_path = "results/$(Dates.format(now(), "Y-mm-ddTHH-MM-SS")).csv"
    end

    # iter -> percentage of infected users per simulation step
    to_return = Dict{Int,Array{Float64,1}}()

    open(output_path, "a") do fhandle
        write(
            fhandle,
            string(
                "sim_step,Δ,δ,c,per_infected,βd,βₑ,βᵢ,γₑ,γₐ,avg_he_size,avg_degree,avg_direct_contacts,",
                "new_users,moved_users,perc_infected_users,perc_infected_locations\n"
                )
        )

        # for randomization purposes
        for iter = 1:niter
            println("Iter $(iter)")

            h = nothing
            added, moved = 0, 0

            # percentage of infected per simulation step
            per_infected_sim = Array{Float64,1}()

            numnodes = length(keys(node_index_map))
            numhes = length(keys(he_index_map))

            # store which users are present in the given timeframe
            usersepoc = zeros(Int, numnodes)

            # evaluation of an initial vector of infected users
            # if it is not given as input 
            if isnothing(vstatus)
                vstatus = fill(1, numnodes)
                vrand = rand(0:100, 1, numnodes)
                for i = 1:numnodes
                    if per_infected  <= vrand[i]
                        vstatus[i] = 0
                    end
                end
            end

            # Initially, all location are safe
            hestatus = zeros(Int, numhes)

            # Storing the new status of each vertex and hyperedge
            vnextstatus = copy(vstatus)
            henextstatus = copy(hestatus)

            # ! check after immunization
            push!(per_infected_sim, sum(vstatus) / length(vstatus))

            # immunization
            istatus = zeros(Int, numnodes)
            irand = rand(0:100, 1, numnodes)

            ihestatus = zeros(Int, numhes)
            iherand = rand(0:100, 1, numhes)

            nextistatus = copy(istatus)
            nextihestatus = copy(ihestatus)
                    
            # Eventually a random set of nodes will start using the tracking app
            # at a certain time
            if app_sim
                Random.seed!(0)
                random_nodes = shuffle(1:numnodes)
                limit = trunc(Int, length(random_nodes) * αₑ)
                random_nodes = random_nodes[1:limit]
                # Keep track of the infected that a node may have been in contact with
                infected_contacts = Dict{Int,Array{Int,1}}()
                map(node -> infected_contacts[node] = [], random_nodes)
                # Nodes that can't spread the infection but can recover
                nodes_in_quarantine = Array{Int,1}()
            end

            ################
            # SIMULATION
            ################
            for t = 1:length(intervals)
                h, added, moved = generatehg(
                                    h, 
                                    df, 
                                    get(intervals, t, 0).first, 
                                    get(intervals, t, 0).second, 
                                    node_index_map, 
                                    he_index_map, 
                                    usersepoc,
                                    t  
                                )   

                # At a specific time  
                if t == 30
                    # appy an immunization strategy, if any, over a percentage of nodes
                    if !(isnothing(nodes_imm_strategy))
                        nodes_to_immunize = nodes_imm_strategy(h)
                        limit = trunc(Int, length(nodes_to_immunize) * αₑ)
                        map(node -> vnextstatus[node] = 0, nodes_to_immunize[1:limit])
                        map(node -> nextistatus[node] = 1, nodes_to_immunize[1:limit])
                    end
                    # appy an immunization strategy, if any, over a percentage of hes
                    if !(isnothing(hes_imm_strategy))
                        hes_to_immunize = hes_imm_strategy(dual(h))
                        limit = trunc(Int, length(hes_to_immunize) * αₑ)
                        map(he -> henextstatus[he] = 0, hes_to_immunize[1:limit])
                        map(he -> nextihestatus[he] = 1, hes_to_immunize[1:limit])
                    end
                end

                if app_sim && t >= 30
                    # For each app user
                    for node in random_nodes
                        # not in quarantine
                        if !(node in nodes_in_quarantine)
                            # find the infected who have been in the same place 
                            neighbours = []
                            hyperedges = collect(keys(gethyperedges(h, node)))
                            for he in hyperedges
                                neighbours = union(neighbours, collect(keys(getvertices(h, he))))
                            end
                            infected_neighoburs = filter(neighbour -> vstatus[neighbour] == 1, neighbours)
                            infected_contacts[node] = union(infected_contacts[node], infected_neighoburs)
                            # Possibly decide to quarantine
                            if rand(1)[1]  < 1 - ℯ^(-0.5 * length(infected_contacts[node]))
                                push!(nodes_in_quarantine, node)
                                continue
                            end
                        end
                    end
                end

                # Estimation of the parameter c
                # based on the distribution
                # of the hyperedge size
                if isnothing(c)
                    dist = Array{Int,1}()
                    
                    for he = 1:nhe(h)
                        push!(dist, length(getvertices(h, he)))
                    end
                    
                    c = median(dist)
                    println(t, " -- ", c)
                end


                #################################
                # Evaluating some stats for 
                # the current hg
                #################################

                # hyperedges average size 
                avg_he_size = .0
                for he = 1:nhe(h)
                    avg_he_size += length(getvertices(h, he))
                end
                avg_he_size /= nhe(h)

                # nodes average degree
                avg_degree = .0
                for v = 1:nhv(h)
                    avg_degree += length(gethyperedges(h, v))
                end
                avg_degree /= nhv(h)

                # number of infected locations with
                # at least two users
                infected_locations = 0
                for he = 1:nhe(h)
                    if hestatus[he] == 1 && length(getvertices(h, he)) > 1
                        infected_locations += 1
                    end
                end

                ########################
                # DIFFUSION ALGORITHM
                ########################

                #
                # PHASE 1 - Agent-to-Environment
                #
                for he = 1:nhe(h)
                    # If the location is immunized,
                    # it cannot spread the infection anymore
                    ihestatus[he] == 1 && rand(1:10) <= 8 && continue
                    # ihestatus[he] == 1 && continue

                    # If the location has at least two users
                    # and it is not infected, it may become contamined.
                    presents = collect(keys(getvertices(h, he)))
                    # Eventually discard nodes in quarantine
                    if app_sim && t >= 30
                        presents = filter(node -> !(node in nodes_in_quarantine), presents)
                    end
                    if length(presents) > 1 && hestatus[he] == 0
                        i = infected(h, he, vstatus, istatus)
                        if rand(1)[1] <  1 - ℯ^- (βₑ * f(i, c))
                            # the location is contaminated
                            henextstatus[he] = 1
                        end
                    elseif rand(1)[1] <  1 - ℯ^- γₑ
                        # the location has been sanitized
                        henextstatus[he] = 0
                    end
                end

                #
                # PHASE 2 - Agent-to-Agent
                #
                avg_direct_contacts = 0
                for v = 1:nhv(h)
                    # an immunized user can't be infected
                    istatus[v] == 1 && continue
                    
                    # if the user is present in the current timeframe
                    if usersepoc[v] == 1
                        # and eventually is not in quarantine
                        if !(app_sim && t >= 30 && v in nodes_in_quarantine)
                            i = 0
                            for he in gethyperedges(h, v)
                                for u in getvertices(h, he.first)
                                    if v != u.first
                                        # if u and v have been together in the same place
                                        # in a time interval less than δ
                                        # then it counts ad a direct contact
                                        if abs(h[v, he.first] - h[u.first, he.first]) <= δ.value
                                            if vstatus[v] == 0
                                                i += vstatus[u.first]
                                            end
                                            avg_direct_contacts += 1
                                        end
                                    end
                                end
                            end
                            # a user becomes infected according to
                            # the following probability
                            if vstatus[v] == 0 && rand(1)[1] < 1 - ℯ^- (βd * i)
                                vnextstatus[v] = 1
                            end
                        end
                    end
                end

                avg_direct_contacts \= sum(usersepoc)


                #
                # PHASE 3 - Environment-to-Agent
                #
                for v = 1:nhv(h)
                    istatus[v] == 1 && continue

                    # if the user is present in the current timeframe
                    if usersepoc[v] == 1

                        # if the user is healthy
                        # and eventually is not in quarantine
                        if vstatus[v] == 0 && !(app_sim && t >= 30 && v in nodes_in_quarantine)
                                i = 0
                                for he in gethyperedges(h, v)
                                    if length(getvertices(h, he.first)) > 1
                                        if ihestatus[he.first] == 0
                                            i += hestatus[he.first]
                                        end
                                    end
                                end
                                if rand(1)[1] < 1 - ℯ^-(βᵢ * f(i, c))
                                    vnextstatus[v] = 1
                                end
                        elseif rand(1)[1] < 1 - ℯ^- γₐ 
                                # the user spontaneously returns healthy
                                vnextstatus[v] = 0
                        end
                    end
                end

                vstatus = copy(vnextstatus)
                istatus = copy(nextistatus)

                hestatus = copy(henextstatus)
                ihestatus = copy(nextihestatus)

                push!(per_infected_sim, sum(vstatus) / length(vstatus))

                to_write = string(
                    "$(t),$(Δ),$(δ),$(c),$(per_infected),$(βd),$(βₑ),$(βᵢ),$(γₑ),$(γₐ),",
                    "$(avg_he_size),$(avg_degree),$(avg_direct_contacts),$(added),$(moved),",
                    "$(sum(vstatus) / length(vstatus)),$(sum(hestatus) / length(hestatus))\n"
                    )

                write(fhandle, to_write)
            end

            push!(to_return, iter => per_infected_sim)
        end
    end

    to_return
end