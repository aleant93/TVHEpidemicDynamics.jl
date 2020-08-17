"""
    infected(h, he, vstatus, istatus)

Count the number of infected nodes within an hyperedge.
"""
function infected(h::Hypergraph, he::Int, vstatus::Array{Int, 1}, istatus::Array{Int, 1}; quarantine::Union{Array{Int, 1}, Nothing} = nothing)
    vertices = getvertices(h, he)
    sum = 0
    for v in vertices
        # if the node is not immunized
        if istatus[v.first] == 0
            if isnothing(quarantine)
                sum += vstatus[v.first]
            else
                if quarantine[v.first] == 0
                    sum += vstatus[v.first]
                end
            end
        end
    end
    sum
end


"""
    f(num_infected, c)

Non-linear function used to bound the infection pressure
for large values of 𝑛.
"""
function f(num_infected::Int, c::Int)
    num_infected > c ? c : num_infected
end
