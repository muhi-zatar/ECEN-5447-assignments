## Dynamic Line Modeling
- Current into/out of each node treated as a parameter in this case, since we're not modeling the generators here
    - If we model the generator, then current becomes a differentiable variable and generator EMF becomes the new parameter
- Solve power flow to get initial conditions for bus voltage (but have to convert this to DQ reference frame)
- Power flow solution also gives initial conditions for current into/out of each node