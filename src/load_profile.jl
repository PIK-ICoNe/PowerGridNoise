"""
    load_profile_model(tspan; γ::Float64 = 1.0, ϵ::Float64 = 0.1, μ_MB::Float64 = 0.0)

A data driven model for a load profile for Residential Electric Power Consumption introduced in [2].

[2] M. Anvari et al., "Data-Driven Load Profiles and the Dynamics of Residential Electric Power Consumption", Arxiv, 2020

This function generates possible realizations of load curves which have the same statistical properties as measured load curves.
"""
function load_profile_model(tspan; γ::Float64 = 1.0, ϵ::Float64 = 0.1, μ_MB::Float64 = 0.0)
    p = [γ, ϵ]
    
    u0 = zeros(3)

    prob_sde_load = SDEProblem(ornstein_uhlenbeck, ornstein_uhlenbeck_ϵ, u0, tspan, p)
    sol = solve(prob_sde_load)

    P_fluc = sqrt.(sol[1, :].^2 + sol[2, :].^2 + sol[3, :].^2) .+ μ_MB

    return P_fluc, sol.t
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