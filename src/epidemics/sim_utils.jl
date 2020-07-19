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

function qinfected(h, he, vstatus, istatus, quarantine)
    vertices = getvertices(h, he)
    sum = 0
    for v in vertices
        # if the node is not immunized
        if istatus[v.first] == 0 && quarantine[v.first] == 0
            sum += vstatus[v.first]
        end
    end
    sum
end
