using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using Revise
using PowerGridNoise, Plots, LaTeXStrings
default(grid = false, foreground_color_legend = nothing, bar_edges = false,  lw=3, framestyle =:box, msc = :auto, dpi=300, legendfontsize = 11, labelfontsize = 15, tickfontsize = 10)

##
sol = pg_frequency_model((0.0, 10.0))
plot(sol, vars = 2, yaxis = L"\omega [rad/s]", xaxis = L"t[s]", legend = false)