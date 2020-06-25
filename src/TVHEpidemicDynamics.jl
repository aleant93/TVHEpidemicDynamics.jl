module TVHEpidemicDynamics

using CSV
using DataFrames
using Dates
using SimpleHypergraphs
using Random

export evaluateintervals!, evaluatedensity!
export evaluatedistance!, evaluate_direct_contacts_distribution!
export evaluate_location_distribution!

export generate_model_data
export buildhg, generatehg
export TVHSIS
export uniform, contacts_based, centrality, acquaintance

include("utils/loader.jl")
include("utils/builder.jl")
include("utils/helper.jl")

include("epidemics/TVHSIS.jl")
include("epidemics/strategies.jl")
include("experiments/compare_strategies.jl")
include("experiments/test_strategies.jl")

end # module
