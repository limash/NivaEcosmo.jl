# Videos of Moro's hourly surface snapshots, one monthly part file of `EcosmoSurfaceWriter`
# (`moro.jl`) at a time: temperature, salinity, surface currents, oxygen and pH, one MP4 each,
# beside the file or in `--out`.
#
# Usage:
#     julia --project=fjordsim_model fjordsim_model/surface_videos.jl \
#         ~/FjordSim_results/moro/snapshots_surface_<tag>_part1.jld2 [--out DIR]
#
# Land is grey and a NaN in the water is magenta, so a failed cell cannot pass for land. Each video
# keeps one colour range, the 1st to 99th percentile of its wet values over the file.

using Dates: DateTime, Millisecond, format
using FjordSim: fjord_data_root
using JLD2: jldopen
using NCDatasets: NCDataset
using Statistics: quantile
import CairoMakie as M

# `moro.jl`'s `START_DATE`: the files store seconds of model time, not dates.
const START_DATE = DateTime(2023, 1, 1)
const FRAMES_PER_SECOND = 24
# Records whose values set a video's colour range, spread evenly over the file.
const COLOUR_RANGE_RECORDS = 50
# One current arrow every this many cells along each axis, and its length in degrees of latitude
# per m s⁻¹.
const ARROW_STRIDE = 12
const ARROW_LENGTH = 0.08

surface_values(name) = (file, iteration) -> file["timeseries/$name/$iteration"][:, :, 1]

# `u` and `v` sit on the cell faces; their means over each cell's two faces are on its centre.
function centred_velocities(file, iteration)
    u = file["timeseries/u/$iteration"][:, :, 1]
    v = file["timeseries/v/$iteration"][:, :, 1]
    return (u[1:end-1, :] .+ u[2:end, :]) ./ 2, (v[:, 1:end-1] .+ v[:, 2:end]) ./ 2
end

function current_speed(file, iteration)
    u, v = centred_velocities(file, iteration)
    return hypot.(u, v)
end

const VIDEOS = (
    (name = "temperature", label = "Temperature (°C)", colormap = :thermal, values = surface_values("T"), arrows = false),
    (name = "salinity", label = "Salinity (g kg⁻¹)", colormap = :haline, values = surface_values("S"), arrows = false),
    (name = "currents", label = "Current speed (m s⁻¹)", colormap = :speed, values = current_speed, arrows = true),
    (name = "oxygen", label = "Oxygen (mmol O₂ m⁻³)", colormap = :viridis, values = surface_values("oxy"), arrows = false),
    (name = "pH", label = "pH", colormap = :plasma, values = surface_values("pH"), arrows = false),
)

"""
    moro_grid()

Cell-centre longitudes and latitudes of the Moro grid and its land mask, from the prepared
`bathymetry.nc`.
"""
function moro_grid()
    longitude, latitude, depth = NCDataset(joinpath(fjord_data_root("moro"), "bathymetry.nc")) do ds
        ds["lon"][:], ds["lat"][:], ds["h"][:, :]
    end
    return (; longitude, latitude, land = depth .>= 0)
end

"""
    colour_range(file, iterations, land, values)

The 1st and 99th percentiles of `values` over the wet, finite cells of `COLOUR_RANGE_RECORDS`
records of the file.
"""
function colour_range(file, iterations, land, values)
    records = unique(round.(Int, range(1, length(iterations); length = COLOUR_RANGE_RECORDS)))
    samples = Float32[]
    for record in records
        field = values(file, iterations[record])
        append!(samples, filter(isfinite, field[.!land]))
    end
    isempty(samples) && return (0.0f0, 1.0f0)
    low, high = quantile(samples, [0.01, 0.99])
    return low < high ? (low, high) : (low - 1, high + 1)
end

function arrow_layer(grid, file, iteration)
    (; longitude, latitude, land) = grid
    cells = [I for I in CartesianIndices(land)[1:ARROW_STRIDE:end, 1:ARROW_STRIDE:end] if !land[I]]
    positions = [M.Point2f(longitude[I[1]], latitude[I[2]]) for I in cells]
    # East in degrees of longitude, so the arrow points the right way on an axis whose aspect is
    # corrected for the convergence of the meridians.
    directions(u, v) = [ARROW_LENGTH * M.Vec2f(u[I] / cosd(latitude[I[2]]), v[I]) for I in cells]
    return positions, directions, centred_velocities(file, iteration)
end

"""
    record_video(path, file, iterations, times, grid, video)

Write one MP4 of `video` over every record of the file.
"""
function record_video(path, file, iterations, times, grid, video)
    (; longitude, latitude, land) = grid
    frame_title(record) = "Moro surface, " * format(START_DATE + Millisecond(round(Int, 1000 * times[record])), "yyyy-mm-dd HH:MM")

    values = M.Observable(video.values(file, first(iterations)))
    title = M.Observable(frame_title(1))
    colorrange = colour_range(file, iterations, land, video.values)

    figure = M.Figure(size = (1100, 800))
    aspect = (last(longitude) - first(longitude)) * cosd((first(latitude) + last(latitude)) / 2) /
             (last(latitude) - first(latitude))
    axis = M.Axis(figure[1, 1]; title, xlabel = "Longitude (°E)", ylabel = "Latitude (°N)",
                  aspect = M.AxisAspect(aspect))
    # Fixed to the grid, so an arrow at the edge is clipped rather than widening the axis.
    M.limits!(axis, extrema(longitude)..., extrema(latitude)...)
    field = M.heatmap!(axis, longitude, latitude, values; colormap = video.colormap, colorrange,
                       nan_color = :magenta)
    M.heatmap!(axis, longitude, latitude, ifelse.(land, 1.0f0, NaN32); colormap = [:gray70, :gray70],
               colorrange = (0, 2), nan_color = :transparent)
    M.Colorbar(figure[1, 2], field; label = video.label)

    if video.arrows
        positions, directions, (u, v) = arrow_layer(grid, file, first(iterations))
        arrows = M.Observable(directions(u, v))
        M.arrows2d!(axis, positions, arrows; color = :black, shaftwidth = 1, tipwidth = 5, tiplength = 5)
    end

    M.record(figure, path, eachindex(iterations); framerate = FRAMES_PER_SECOND) do record
        values[] = video.values(file, iterations[record])
        title[] = frame_title(record)
        video.arrows && (arrows[] = directions(centred_velocities(file, iterations[record])...))
    end

    println("Wrote $path")
    return path
end

function main(args)
    value(flag, default) = (i = findfirst(==(flag), args); isnothing(i) ? default : args[i+1])
    path = first(args)
    out = value("--out", dirname(abspath(path)))
    mkpath(out)

    grid = moro_grid()
    jldopen(path) do file
        iterations = sort(parse.(Int, keys(file["timeseries/t"])))
        times = [file["timeseries/t/$iteration"] for iteration in iterations]
        stem = first(splitext(basename(path)))
        for video in VIDEOS
            record_video(joinpath(out, "$(stem)_$(video.name).mp4"), file, iterations, times, grid, video)
        end
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
