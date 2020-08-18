abstract type Abstract_Simulation_Model end

struct SIS <: Abstract_Simulation_Model end
struct SIS_infected_loc <: Abstract_Simulation_Model end
struct SIS_tracing <: Abstract_Simulation_Model end
struct SIS_sanification <: Abstract_Simulation_Model end
