"""

"""
function simulate_immuni(
        df::DataFrame,
        intervals::Dict{Int, Pair{DateTime, DateTime}},
        user2vertex::Dict{String, Int},
        loc2he::Dict{String, Int},
        δ::Dates.Millisecond;
        Δ::Union{Int,TimePeriod,Nothing} = nothing,
        vstatus::Union{Array{Int, 1}, Nothing} = nothing,
        per_infected::Int = 20,
        c::Union{Int, Nothing} = 5,
        βd::Float64 = 0.2,
        βₑ::Float64 = 0.06,
        βᵢ::Float64 = 0.1,
        γₑ::Float64 = 0.06,
        γₐ::Float64 = 0.1,
        αᵥ::Float64 = 0.0,
        αₑ::Float64 = 0.0,
        αᵢ::Float64 = 0.0,
        lockdown::Bool = false,
        βₗ::Float64 = 0.0,
        imm_start::Int = 0,
        app_strategy::Union{Function, Nothing} = nothing,
        nodes_imm_strategy::Union{Function, Nothing} = nothing,
        hes_imm_strategy::Union{Function, Nothing} = nothing,
        nodes_kwargs::Dict = Dict{}(),
        hes_kwargs::Dict = Dict{}(),
        niter::Int = 1,
        output_path::Union{AbstractString, Nothing} = nothing
)

    if isnothing(output_path)
        if !isdir("results/")
            Filesystem.mkdir("results/")
        end

        output_path = "results/$(Dates.format(now(), "Y-mm-ddTHH-MM-SS")).csv"
    end

    # iter -> percentage of infected users per simulation step
    to_return = Dict{Int, Array{Float64, 1}}()

    # evaluation of an initial vector of infected users
    # if it is not given as input
    if isnothing(vstatus)
        vstatus = fill(1, length(user2vertex))
        vrand = rand(0:100, 1, length(user2vertex))
        for i=1:length(user2vertex)
            if per_infected  <= vrand[i]
                vstatus[i] = 0
            end
        end
    end


    open(output_path, "a") do fhandle
        write(
            fhandle,
            string(
                "sim_step,Δ,δ,c,per_infected,βd,βₑ,βᵢ,γₑ,γₐ,avg_he_size,",
                "avg_degree,avg_direct_contacts,new_users,",
                "moved_users,perc_infected_users,perc_infected_locations\n"
                )
        )
        # for randomization purposes
        for iter=1:niter
            println("Iter $(iter)")

            h = nothing
            added, moved = 0, 0

            # percentage of infected per simulation step
            per_infected_sim = Array{Float64, 1}()

            # store which users are present in the given timeframe
            usersepoc = zeros(Int, length(user2vertex))

            # Initialize the status vector of the nodes
            _vstatus = copy(vstatus)
            # Initially, all location are safe
            hestatus = zeros(Int, length(loc2he))

            # Storing the new status of each vertex and hyperedge
            vnextstatus = copy(vstatus)
            henextstatus = copy(hestatus)

            push!(per_infected_sim, sum(_vstatus) / length(_vstatus))

            ################
            # IMMUNIZATION
            ################
            istatus = rand(0:0, 1, length(user2vertex))
            ihestatus = rand(0:0, 1, length(loc2he))

            nextistatus = copy(istatus)
            nextihestatus = copy(ihestatus)

            ################
            # TRACING
            ################
            to_trace = Dict{Int, Set{Int}}()
            users_to_trace = Set{Int}()
            quarantine = fill(0, length(user2vertex))

            ################
            # SIMULATION
            ################
            for t=1:length(intervals)
                h, added, moved = generatehg!(
                                    h,
                                    df,
                                    get(intervals, t, 0).first,
                                    get(intervals, t, 0).second,
                                    user2vertex,
                                    loc2he,
                                    usersepoc,
                                    t
                                )

                isnothing(h) && continue

                # Estimation of the parameter c
                # based on the distribution
                # of the hyperedge size
                if isnothing(c)
                    dist = Array{Int, 1}()

                    for he=1:nhe(h)
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
                for he=1:nhe(h)
                    avg_he_size += length(getvertices(h, he))
                end
                avg_he_size /= nhe(h)

                # nodes average degree
                avg_degree = .0
                for v=1:nhv(h)
                    avg_degree += length(gethyperedges(h, v))
                end
                avg_degree /= nhv(h)

                # number of infected locations with
                # at least two users
                infected_locations = 0
                for he=1:nhe(h)
                    if hestatus[he] == 1 && length(getvertices(h, he)) > 1
                        infected_locations += 1
                    end
                end


                ########################
                # IMMUNIZATION
                # immunizing a node does not have effect on its health status:
                # a node that has been immunized cannot get sick,
                # but it can still spread the contagion
                ########################

                # start the immunization process
                # in a given timeframe
                if t == imm_start
                    # apply the given immunization strategy
                    # over αᵥ nodes
                    if !isnothing(nodes_imm_strategy)
                        to_immunize = nodes_imm_strategy(h, αᵥ; nodes_kwargs...)
                        map(v -> nextistatus[v] = 1, to_immunize)
                    end

                    # apply the given immunization strategy
                    # over αₑ hyperedges
                    if !isnothing(hes_imm_strategy)
                        to_immunize = hes_imm_strategy(dual(h), αₑ; hes_kwargs...)
                        map(he -> nextihestatus[he] = 1, to_immunize)
                    end

                    if !isnothing(app_strategy)
                        u2trace = app_strategy(h, αᵥ; nodes_kwargs...)
                        map(user -> to_trace[user] = Set{Int}(), u2trace)
                        users_to_trace = keys(to_trace)
                    end
                end


                ########################
                # TRACING
                ########################
                for v in users_to_trace
                    if quarantine[v] == 0
                        i = 0
                        for u in to_trace[v]
                            i += _vstatus[u]
                        end
                        if rand() < 1 - ℯ ^ - (1 * i)
                            quarantine[v] = 1
                        end
                    end
                end

                println("Nodes in quarantine $(sum(quarantine)) at iter $(t)")

                ########################
                # DIFFUSION ALGORITHM
                ########################

                #
                # PHASE 1 - Agent-to-Environment
                #
                for he=1:nhe(h)
                    # If the location is immunized
                    # and the lockdown is active,
                    # it cannot spread the infection anymore
                    ihestatus[he] == 1 && lockdown && continue

                    # If the location has at least two users
                    # and it is not infected, it may become contamined.
                    if length(getvertices(h, he)) > 1 && hestatus[he] == 0
                        i = qinfected(h, he, _vstatus, istatus, quarantine)
                        if rand() <  1 - ℯ ^ - (βₑ * f(i, c))
                            # the location is contaminated
                            henextstatus[he] = 1
                        end
                    elseif rand() <  1 - ℯ ^ - γₑ
                        # the location has been sanitized
                        henextstatus[he] = 0
                    end
                end

                #
                # PHASE 2 - Agent-to-Agent
                #
                avg_direct_contacts = 0
                for v=1:nhv(h)
                    # a node that has been immunized cannot become infected
                    (istatus[v] == 1 || quarantine[v] == 1) && continue

                    # if the user is present in the current timeframe
                    if usersepoc[v] == 1
                        i = 0
                        for he in gethyperedges(h, v)

                            # if the lockdown is active,
                            # a direct contact in that location cannot happen
                            ihestatus[he.first] == 1 && lockdown  && continue

                            for u in getvertices(h, he.first)
                                if v != u.first
                                    # if v and u are using our tracing app
                                    # then we'll save this contact
                                    if v in users_to_trace && u.first in users_to_trace
                                        push!(
                                            get!(to_trace, v, Set{Int}()),
                                            u.first
                                        )

                                        push!(
                                            get!(to_trace, u.first, Set{Int}()),
                                            v
                                        )
                                    end

                                    # if u and v have been together in the same place
                                    # in a time interval less than δ
                                    # then it counts ad a direct contact
                                    if abs(h[v, he.first] - h[u.first, he.first]) <= δ.value
                                        if _vstatus[v] == 0
                                            i += _vstatus[u.first]
                                        end
                                        avg_direct_contacts += 1


                                    end
                                end
                            end
                        end
                        # a user becomes infected according to
                        # the following probbaility
                        if _vstatus[v] == 0 && rand() < 1 - ℯ ^ - (βd * i)
                            vnextstatus[v] = 1
                        end
                    end
                end

                avg_direct_contacts \= sum(usersepoc)


                #
                # PHASE 3 - Environment-to-Agent
                #
                for v=1:nhv(h)

                    # if the user is present in the current timeframe
                    if usersepoc[v] == 1

                        # if the user is healthy and
                        # it is not immunized,
                        # it may become ill
                        if _vstatus[v] == 0 && istatus[v] == 0 && quarantine[v] == 0
                            i = 0
                            i_immunized = 0

                            for he in gethyperedges(h, v)

                                # if the lockdown is active,
                                # an indirect contact in that location
                                # cannot take place
                                ihestatus[he.first] == 1 && lockdown && continue

                                if length(getvertices(h, he.first)) > 1
                                    if ihestatus[he.first] == 0
                                        i += hestatus[he.first]
                                    else
                                        i_immunized += hestatus[he.first]
                                    end
                                end
                            end

                            if rand() < 1 - ℯ ^ -( (βᵢ * i) + ( (βᵢ * (1 - βₗ)) * i_immunized) )
                                vnextstatus[v] = 1
                            end

                        elseif rand() < 1 - ℯ ^ - γₐ
                            # the user spontaneously returns healthy
                            vnextstatus[v] = 0
                        end
                    end
                end

                _vstatus = copy(vnextstatus)
                istatus = copy(nextistatus)

                hestatus = copy(henextstatus)
                ihestatus = copy(nextihestatus)

                #push!(per_infected_sim, sum(_vstatus) / (length(_vstatus) - sum(istatus)))
                push!(per_infected_sim, sum(_vstatus) / length(_vstatus))

                to_write = string(
                    "$(t),$(Δ),$(δ),$(c),$(per_infected),$(βd),$(βₑ),$(βᵢ),$(γₑ),$(γₐ),",
                    "$(avg_he_size),$(avg_degree),$(avg_direct_contacts),$(added),$(moved),",
                    "$(sum(_vstatus)/length(_vstatus)),$(sum(hestatus)/length(hestatus))\n"
                    )

                write(fhandle, to_write)
            end

            push!(to_return, iter=>per_infected_sim)
        end
    end

    to_return
end
