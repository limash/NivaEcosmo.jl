# Compare the Julia box run with the Fortran (pyfabm) one, variable by variable.
#
# Usage:
#     julia --project=box_model box_model/compare.jl <fortran output.nc> [<julia output.nc>]
#
# For every variable the two files share, prints the largest difference over the year relative to
# that variable's largest magnitude, worst first, and writes compare.png next to the Julia output
# with the two runs overlaid for the main state variables.

using NCDatasets
import CairoMakie as M

const PLOTTED = ("ECO_no3", "ECO_nh4", "ECO_pho", "ECO_sil", "ECO_oxy", "ECO_fla", "ECO_dia",
                 "ECO_flachl", "ECO_diachl", "ECO_microzoo", "ECO_mesozoo", "ECO_det",
                 "ECO_dom", "ECO_cal", "ECO_sed1", "CO2_dic", "CO2_alk", "CO2_pH")

read_variables(path) = NCDataset(path) do ds
    Dict(name => Float64.(ds[name][:]) for name in keys(ds) if name != "time")
end

"""
    relative_difference(julia, fortran)

The largest absolute difference over time, scaled by the largest magnitude either run reaches.
"""
function relative_difference(julia, fortran)
    scale = max(maximum(abs, julia), maximum(abs, fortran))
    return scale == 0 ? 0.0 : maximum(abs.(julia .- fortran)) / scale
end

function main(args)
    isempty(args) && error("usage: compare.jl <fortran output.nc> [<julia output.nc>]")
    fortran_path = args[1]
    julia_path = length(args) > 1 ? args[2] : joinpath(@__DIR__, "output.nc")

    fortran = read_variables(fortran_path)
    julia = read_variables(julia_path)

    shared = sort(collect(intersect(keys(fortran), keys(julia))))
    missing_names = setdiff(keys(fortran), keys(julia))
    differences = sort([(name, relative_difference(julia[name], fortran[name])) for name in shared]; by = last, rev = true)

    println("Largest difference relative to the variable's largest magnitude (Julia vs Fortran):")
    for (name, difference) in differences
        println(rpad(name, 60), difference)
    end
    isempty(missing_names) || println("\nOnly in the Fortran output: ", join(sort(collect(missing_names)), ", "))

    fig = M.Figure(size = (1400, 1400))
    for (n, name) in enumerate(PLOTTED)
        ax = M.Axis(fig[cld(n, 3), mod1(n, 3)]; title = name, xlabel = "day")
        days = 0:length(fortran[name])-1
        M.lines!(ax, days, fortran[name]; label = "Fortran", linewidth = 3, color = :gray70)
        M.lines!(ax, days, julia[name]; label = "Julia", linestyle = :dash, color = :black)
        n == 1 && M.axislegend(ax; labelsize = 10)
    end
    path = joinpath(dirname(julia_path), "compare.png")
    M.save(path, fig)
    println("\nWrote $path")

    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
