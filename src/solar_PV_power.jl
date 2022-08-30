"""
    PV_power_model(tspan; α::Float64 = 1.0,  D_diff_0::Float64 = 0.1, σ_ϵ2::Float64 = 0.5, λ::Float64 = 1.0, u0)

Generates an intermittent wind power fluctuation time series by use of the Langevin-type model introduced in [3].

## Arguments:
- `tspan::NTuple{2, Float64}`: Time span of the power time series.
    
## Keywords:
- `α::Float64`: Steepness of the 
- `D_diff::Float64`: Diffusion Term / Second KM coefficient
- `σ_ϵ2::Float64:` Jump Amplitude
- `λ::Float64`: Jump rate
- `u0::Float64`: Initial condition
"""
function solar_PV_power_model(tspan, u0; cloudy::Bool, sunny::Bool, α::Float64 = 0.021,  D_diff_0::Float64 = 0.15, symmetry::Bool = false, σ_ϵ2::Float64 = 0.5, λ::Float64 = 1.0)
    p = [α, D_diff_0, σ_ϵ2, λ, symmetry, cloudy, sunny] # Create the parameter array

    prob_sde_sol = SDEProblem(clear_sky_deterministic, clear_sky_noise, u0, tspan, p)
    sol = solve(prob_sde_sol)

    return sol
end

function clear_sky_deterministic(du, u, p, t)
    α, _, _, _, symmetry, cloudy, sunny = p

    du[1] = u[2] + u[3]
    du[2] = D_drift(x = u[1],α = α, symmetry = symmetry, cloudy = cloudy, sunny = sunny)
    du[3] = 0.0 # Diffusion
    return nothing
end

function clear_sky_noise(du, u, p, t)    
    D_diff_0 = p[2]
    du[1] = 0.0 
    du[2] = 0.0 # Drift
    du[3] = sqrt(D_diff(x = u[1], D_diff_0 = D_diff_0)) 
    return nothing
end

"""
    D_drift(x, α, symmetry, cloudy, sunny)

Function for the Drift Term / First Kramers-Moyal (KM) coefficient.
"""
function D_drift(;x, α::Float64, symmetry, cloudy, sunny)
    D_drift = α * (symmetry * (-32 * x^3 + 4 * x) + (symmetry - 1) * (-40 * x^3 + 4 * x + 0.3) + cloudy * (-20 * (x + 0.2)) + sunny * (-20 * (x - 0.2)))
    return D_drift
end

"""
    D_diff(x, D_diff_0)

Function for the Diffusion Term / Second Kramers-Moyal (KM) coefficient.
"""
function D_diff(;x, D_diff_0)
    D_drift = D_diff_0 * sqrt(abs(x * (x + 0.4)))

    return D_drift
end