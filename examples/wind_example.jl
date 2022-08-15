using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using Revise
using Random
using PowerGridNoise, Plots, LaTeXStrings
default(grid = false, foreground_color_legend = nothing, bar_edges = false,  lw=3, framestyle =:box, msc = :auto, dpi=300, legendfontsize = 11, labelfontsize = 15, tickfontsize = 10)

##
D = 0.1

tspan = (0.0, 200.0)

sol = wind_power_model(tspan, D = D)

plot(sol, idxs = 1, legend = false, xlabel = L"t[s]", ylabel = L"x(t)")