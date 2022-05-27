"""
    pg_frequency_model(tspan, P0::Float64 = 0.001641, ϵ::Float64 = 0.00105, c1::Float64 = 0.008311, c2::Float64 = 0.000030)

A Data-Driven Model of the Power-Grid Frequency Dynamics which was introduced in [1].

[1] L. R. Gorjão et al., "Data-Driven Model of the Power-Grid Frequency Dynamics," in IEEE Access, vol. 8, pp. 43082-43097, 2020, doi: 10.1109/ACCESS.2020.2967834.

This function generates a random frequency trajectory which features the same statistical properties as real, measured power grid frequency time series.
The default parameters [c1, c2, P0] are taken from [1] and depict the central european case (CE).

## Arguments:
- `tspan::NTuple{2, Float64}`: Time span of the frequency time series.

## Keywords:
- `P0::Float64`: Initial Power Imbalance
- `ϵ::Float64`: White Noise Amplitude
- `c1::Float64`: Primary Control Parameter
- `c2::Float64`: Secondary Control Parameter
"""
function pg_frequency_model(tspan::NTuple{2, Float64}; P0::Float64 = 0.001641, ϵ::Float64 = 0.00105, c1::Float64 = 0.008311, c2::Float64 = 0.000030)
    new_dispatch = 60 * 15 # A new dispatch occurs every 15 minutes
    n_dispatch = floor(Int, tspan[end] / new_dispatch) # Number of dispatches in tspan
    dispatch_times = collect(range(0, step = new_dispatch, length = n_dispatch + 1)) 
    p = vcat([c1, c2, P0, ϵ], dispatch_times) # Simulation parameters

    dt = 0.01 # fixed step size

    #pg_ω_deter_step = pg_ω_deter_step(du, u, p, t)            # Deterministic part of the frequency dynamics: steps in power due to new dispatches
    #pg_ω_noise = pg_ω_noise(du,u,p,t)                         # Stochastic part of the frequency dynamics
    #dispatch_condition = dispatch_condition(u, t, integrator) # Dispatches only occur during the scheduled times

    dispatch = DiscreteCallback(dispatch_condition, affect!)  # Callback which checks if a dispatch should occur and samples a new power mismatch at the dispatch
    prob_sde = SDEProblem(pg_ω_deter_step, pg_ω_noise, [0.0, 0.0], tspan, p, Noise = WienerProcess(0.0, 0.0, 0.0))
    sol = solve(prob_sde, tstops = dispatch_times, callback = dispatch, dt = dt, adaptive = false)

    return sol #, t -> (sol(t)[2], sol(t, Val{1})[2]) # ω and dω
end

"""
    pg_ω_deter_step(du, u, p, t)

Deterministic part of the frequency dynamics. The power jumps due to the new dispatches with occur every 15 minutes. 
This function is only used internally during the simulation of the stochastic differential equation.
"""
function pg_ω_deter_step(du, u, p, t)
    θ, ω = u[1], u[2]
    dθ, dω = du[1], du[2]
    
    # ΔP jumps at every new dispatch
    c1, c2, ΔP = p # New simulation parameters
    dθ = ω
    dω = -c1 * ω -c2 * θ + ΔP
end

"""
    pg_ω_noise(du,u,p,t)

Stochastic part of the frequency dynamics.
This function is only used internally during the simulation of the stochastic differential equation.
"""
function pg_ω_noise(du, u, p, t)
    ϵ = p[4]

    du[1] = 0.0
    du[2] = ϵ # Only the ω dynamics has an additional white noise term
end

"""
    dispatch_condition(u, t, integrator)

Function which checks if a new dispatch should occur.
This function is only used internally during the simulation of the stochastic differential equation.
"""
function dispatch_condition(u, t, integrator)
    dispatch_times = integrator.p[end]
    t ∈ dispatch_times
end

"""
    affect!(integrator)

Function which samples the new power mismatch after each dispatch.
This function is only used internally during the simulation of the stochastic differential equation.
"""
function affect!(integrator)
    # integrator.u[3] = rand(Uniform(0.0001, 0.002)) 
    integrator.p[3] = rand(Uniform(-1 * 10^-3, 1 * 10^-3))  # New Power Imbalance
end