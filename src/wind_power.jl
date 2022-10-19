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
- `fluc_bound::Float64`: Limit which the fluctuation should not exceed. Everything above this is cut off.
- `Δt::Float64`: Swing in time of the dynamical process
"""
function wind_power_model(tspan::NTuple{2, Float64}; D::Float64 = 0.1, γ::Float64 = 1.0, g::Float64 = 0.5, x0::Float64 = 2.0, ϵ::Float64 = 1.0, fluc_bound::Float64 = 1.0, Δt::Float64 = 1000.0)
    p = [g, x0, γ, D, ϵ] # Create the parameter array
    u0 = [g * x0, 0.0] # Initial condition
    tspan = (tspan[1], tspan[2] + Δt)

    prob_sde_wind = SDEProblem(wind_model_deterministic, wind_model_noise, u0, tspan, p, noise = WienerProcess(0.0, 0.0, 0.0), tstops = Δt)
    sol = solve(prob_sde_wind)

    swing_in_t_idx = findfirst(sol.t .>= Δt) 

    mean_x = g * x0
    x = sol[1, :] .- mean_x # Shifting the mean of the time series to 0.0,

    # Cutting of any fluctuation that becomes bigger than the allowed fluctuation bound
    upper_bound_idx = findall(x .> fluc_bound) 
    x[upper_bound_idx] .= fluc_bound

    lower_bound_idx = findall(x .< -fluc_bound)
    x[lower_bound_idx] .= -fluc_bound

    # Throw away the part of the time series before Δt
    x = x[swing_in_t_idx:end]
    t = sol.t[swing_in_t_idx:end] .- Δt

    return x, t
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