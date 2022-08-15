"""
    wind_power_model(D::Float64, γ::Float64 = 1.0, g::Float64 = 0.5, x0::Float64 = 2.0)

Generates an intermittent wind power fluctuation time series by use of the Langevin-type model introduced in [3].

[3] K. Schmietendorf et al., "The impact of turbulent renewable energy production on power grid stability and quality", EUROPEAN PHYSICAL JOURNAL B, 2017

## Arguments:
- `tspan::NTuple{2, Float64}`: Time span of the power time series.
    
## Keywords:
- `D::Float64`: Intermittence strength
- `γ::Float64`: Damping coefficient
- `g::Float64`
- `x0::Float64`
"""

function wind_power_model(tspan; D::Float64, γ::Float64 = 1.0, g::Float64 = 0.5, x0::Float64 = 2.0, ϵ::Float64 = 1.0)
    p = [g, x0, γ, D, ϵ] # Create the parameter array
    u0 = zeros(2)

    prob_sde_wind = SDEProblem(wind_model_deterministic, wind_model_noise, u0, tspan, p, noise = WienerProcess(0.0, 0.0, 0.0))
    sol = solve(prob_sde_wind)

    return sol
end

function wind_model_deterministic(du, u, p, t)
    g, x0, γ, D, ϵ = p

    du[1] = u[1] * (g - (u[1] / x0)) + sqrt(D * u[1]^2) * u[2]
    du[2] = -γ * u[2]

    return nothing
end

function wind_model_noise(du, u, p, t)
    ϵ = p[5]
    
    du[1] = 0.0 
    du[2] = ϵ

    return nothing
end