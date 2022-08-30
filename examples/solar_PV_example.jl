using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using Revise
using Random
using PowerGridNoise, Plots, LaTeXStrings
default(grid = false, foreground_color_legend = nothing, bar_edges = false,  lw=3, framestyle =:box, msc = :auto, dpi=300, legendfontsize = 11, labelfontsize = 15, tickfontsize = 10)

##
tspan = (0.0, 100.0)
u0 = zeros(3)

sol = solar_PV_power_model(tspan, u0; cloudy = false, sunny = false)

plot(sol, idxs = 1)

##
import PowerGridNoise.D_drift
Z = collect(0.02:0.001:1.0)

D_d = map(x -> D_drift(x = x, α = 0.021, symmetry = 0.0, cloudy = 0.0, sunny = 0.0), Z)

plot(Z, D_d)

##
import PowerGridNoise.D_diff

D_dif = map(x -> D_diff(x = x, D_diff_0 = 0.15), Z)

plot(Z, D_dif)