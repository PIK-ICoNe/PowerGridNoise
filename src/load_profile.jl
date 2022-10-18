"""
    load_profile_model(tspan; γ::Float64 = 1.0, ϵ::Float64 = 0.1, μ_MB::Float64 = 0.0)

A data driven model for a load profile for Residential Electric Power Consumption introduced in [2].

[2] M. Anvari et al., "Data-Driven Load Profiles and the Dynamics of Residential Electric Power Consumption", Arxiv, 2020

This function generates possible realizations of load curves which have the same statistical properties as measured load curves.
The default parameters [start_demand, ϵ, γ, μ_MB] are taken from [2] and and were extracted from the NOVAREF data from 2018.

## Arguments:
- `tspan::NTuple{2, Float64}`: Time span of the frequency time series.
    
## Keywords:
- `start_demand::Float64`: Initial demand
- `ϵ::Float64`: Noise Amplitude for Ornstein Uhlenbeck
- `γ::Float64`: Damping coefficient for Ornstein Uhlenbeck
- `μ_MB::Float64`: Observed shift from zero
- `Δt::Float64`: Swing in time of the dynamical process
"""
function load_profile_model(tspan::NTuple{2, Float64}; start_demand = 42.89173333414218, γ::Float64 = 0.015581760061431082, ϵ::Float64 = 33.806740449173745, μ_MB::Float64 = 0.026988400097746185, Δt::Float64 = 1000.0)
    p = [γ, ϵ] # create the parameter array
    tspan = (tspan[1], tspan[2] + Δt)

    u0 = start_demand * ones(3) # sets the initial condition for each OU process

    prob_sde_load = SDEProblem(ornstein_uhlenbeck, ornstein_uhlenbeck_ϵ, u0, tspan, p, tstops = Δt)
    sol = solve(prob_sde_load)

    P_fluc = sqrt.(sol[1, :].^2 + sol[2, :].^2 + sol[3, :].^2) .+ μ_MB # Equation 3 in paper [2]

    swing_in_t_idx = findfirst(sol.t .>= Δt) 
    
    # Throw away the part of the time series before Δt
    P_fluc = P_fluc[swing_in_t_idx:end] 
    t = sol.t[swing_in_t_idx:end] .- Δt

    return P_fluc, t
end

"""
    ornstein_uhlenbeck(du, u, p, t)

Defines the deterministic part of the Ornstein Uhlenbeck process. Sources [2] states that three processes capture the load dynamics best, hence we also use three.
This function is only used internally during the simulation of the stochastic differential equation.
"""
function ornstein_uhlenbeck(du, u, p, t)
    γ = p[1]

    du[1] = -γ * u[1]
    du[2] = -γ * u[2]
    du[3] = -γ * u[3]
end

"""
    ornstein_uhlenbeck_ϵ(du, u, p, t)

Defines the stochastic part of the Ornstein Uhlenbeck process. ϵ defines the noise amplitude.
This function is only used internally during the simulation of the stochastic differential equation.
"""
function ornstein_uhlenbeck_ϵ(du, u, p, t)
    ϵ = p[2]

    du[1] = ϵ
    du[2] = ϵ
    du[3] = ϵ
end