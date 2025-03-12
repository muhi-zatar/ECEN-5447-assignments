```
This is a network model module, currently not in use, and will not be used.
This is as a placeholder for the line network when ready.
```

# Defining the module
module NetworkModel

# Exporing the functions
export create_three_bus_network

using ..PowerFlowSolver

# Function to create a 3-bus network with G1 as infinite source, G2 as our detailed model, and a load
function create_three_bus_network(;
    base_MVA = 100.0,
    # Generator 2 parameters
    G2_P = 1.0,    # pu
    G2_V = 1.05,   # pu
    # Load parameters
    load_P = 1.5,  # pu
    load_Q = 0.5   # pu
)
    # Create buses
    buses = [
        # Bus 1: Slack bus (infinite source G1)
        Bus(1, SLACK, 1.06, 0.0, 0.0, 0.0, 0.0, 0.0),
        
        # Bus 2: PV bus (generator G2 with our detailed model)
        Bus(2, PV, G2_V, 0.0, G2_P, 0.0, 0.0, 0.0),
        
        # Bus 3: PQ bus (load)
        Bus(3, PQ, 1.0, 0.0, 0.0, 0.0, load_P, load_Q)
    ]
    
    # Create branches (transmission lines)
    branches = [
        # Line from Bus 1 to Bus 2
        Branch(1, 2, 0.01, 0.1, 0.02),
        
        # Line from Bus 2 to Bus 3
        Branch(2, 3, 0.02, 0.2, 0.04),
        
        # Line from Bus 3 to Bus 1
        Branch(3, 1, 0.03, 0.3, 0.06)
    ]
    
    # Create and return the network
    return Network(buses, branches, base_MVA)
end

end # module