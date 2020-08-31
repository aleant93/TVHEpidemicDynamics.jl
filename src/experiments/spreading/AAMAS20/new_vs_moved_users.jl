using CSV
using PyPlot

"""
 Counts how many new users come into play at each simulation step,
 along with how many users change their location.
"""

# each log file is ok to evaluate this information
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
