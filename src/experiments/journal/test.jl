using SimpleHypergraphs
using TVHEpidemicDynamics

h = Hypergraph{Int}(6,4)

h[1:3, 1] .= 1
h[2:4, 2] .= 2
h[2:5, 3] .= 3
h[1, 4] = 4
h[3, 4] = 4
h[6, 4] = 4

d = [length(gethyperedges(h, v)) for v in 1:nhv(h)]

to_return = sort!(collect(1:nhv(h)), by = x -> d[x], rev = true)
to_return[1:ceil(Int, length(to_return)*α)]

R, P = hrwr(h; heweights=nothing, α =.8, ϵ = .0001, iter=1000)

rank = R[1, :]

sort!(collect(1:nhv(h)), by = x -> rank[x], rev = true)

data = [(v, R[1, v]) for v in 1:nhv(h)]
sort!(data, by = x -> x[2], rev = true)
to_return = [pair[1] for pair in data]


using CSV
using PyPlot

df = CSV.read("src/experiments/journal/results/csv/immuni415_V-tracing-60_2020-07-17T10-54-17.csv")

clf()
figure(figsize=(7,4))

plot(df[!, :new_users])
plot(df[!, :moved_users])

legend(["New users", "Moved users"], fontsize="large", ncol=2)

xlabel("Time intervals", fontweight="semibold", labelpad=10, fontsize="x-large")
ylabel("Number of users", fontweight="semibold", fontsize="x-large", labelpad=10)

gcf()

plt.tight_layout(.5)
savefig("mv_vs_new.png")