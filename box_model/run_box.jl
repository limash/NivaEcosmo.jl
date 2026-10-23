# Single-box (0D) model of Drammensfjorden with the Julia ECOSMO: the twin of niva-ecosmo's
# `box_model/run_box.py`, which runs the Fortran model through pyfabm.
#
# The box is well mixed, with prescribed temperature, salinity, light and wind (forcing.csv).
# Nutrients are supplied by exchange with the river and the outer fjord (boundary.yaml):
#
#     dC/dt = ECOSMO rates (interior + surface + bottom) + (Q_river/V) (C_river - C) + (Q_ocean/V) (C_ocean - C)
#
# Everything the Python driver does is done here the same way and in the same order — forward
# Euler, FABM's clipping of each state variable to its minimum after the step, the Knudsen
# exchange, the diel light cycle, numpy's linear interpolation of the daily forcing — so that the
# two runs differ only by how the two languages round, and `compare.jl` can say by how much.
#
# Usage:
#     julia --project=box_model box_model/run_box.jl             # full run -> output.nc, output.png
#     julia --project=box_model box_model/run_box.jl --closed    # no exchange, no burial: P conservation check

using NivaEcosmo
using NivaEcosmo: interior_rates, bottom_fluxes, oxygen_surface_flux, carbonate_diagnostics,
    carbon_dioxide_surface_flux, day_length, gotm_light, attenuation_coefficient, ECOSMO_MINIMA,
    ECOSMO_PARAMETER_NAMES
using Dates
using NCDatasets
using YAML
import CairoMakie as M

const HERE = @__DIR__
# The inputs, shared with the tests and `fjordsim_model/`
const DATA = joinpath(HERE, "..", "data")

# Box geometry and constant forcing, as in run_box.py
const LATITUDE = 59.7          # degrees N
const DEPTH = 10.0             # box thickness H (m)
const AREA = 50.0e6            # box surface area (m2)
const VOLUME = AREA * DEPTH    # (m3)
const S_OCEAN = 30.0           # salinity of inflowing outer-fjord water (PSU), for the Knudsen exchange
const Q_OCEAN = nothing        # constant outer-fjord inflow (m3/s); nothing = derive from the salt balance
const BOTTOM_STRESS = 0.0      # (Pa) below crBotStr, so that detritus settles to the sediment
const PCO2_AIR = 420.0         # atmospheric CO2 mole fraction (ppm)
const DT = 600.0               # time step (s)
const DIEL_LIGHT = true        # impose a diel cycle on the daily-mean shortwave radiation

# boundary.yaml units (µmol/L) to ECOSMO units (mgC/m3)
const TO_MODEL_UNITS = Dict("no3" => 6.625 * 12.01, "nh4" => 6.625 * 12.01, "pho" => 106.0 * 12.01, "sil" => 6.625 * 12.01)

const INTERIOR = keys(ECOSMO_MINIMA)
const SEDIMENT = (:sed1, :sed2, :sed3, :sed4)

# Phosphorus-bearing pools (all in mgC Redfield units), for the --closed conservation check
const P_INTERIOR = (:pho, :fla, :dia, :microzoo, :mesozoo, :det, :dom)
const P_BOTTOM = (:sed1, :sed3)
const BURIAL_PARAMETERS = (:burialRt, :burialRt2, :burialRtO, :burialRt2O, :burialRtC, :burialRt2C)

# Names in the output file, as pyfabm writes them
output_name(name) = name in (:dic, :alk) ? "CO2_$name" : "ECO_$name"

#####
##### Forcing
#####

"""
    read_forcing(path)

The daily forcing table: a `NamedTuple` of the `date` column and one `Float64` vector per other
column. Lines starting with `#` are comments.
"""
function read_forcing(path)
    lines = filter(line -> !isempty(line) && !startswith(line, "#"), readlines(path))
    header = Symbol.(split(lines[1], ","))
    rows = [split(line, ",") for line in lines[2:end]]

    columns = map(enumerate(header)) do (n, name)
        name == :date ? [Date(row[n]) for row in rows] : [parse(Float64, row[n]) for row in rows]
    end

    return NamedTuple{Tuple(header)}(Tuple(columns))
end

"""
    numpy_interp(x, xp, fp)

`numpy.interp` for a scalar `x`, with its exact arithmetic: clamped at both ends, the left node's
value at a node, and `slope * (x - xp[j]) + fp[j]` in between.
"""
function numpy_interp(x, xp, fp)
    x < xp[1] && return fp[1]
    x >= xp[end] && return fp[end]
    j = searchsortedlast(xp, x)
    xp[j] == x && return fp[j]
    slope = (fp[j+1] - fp[j]) / (xp[j+1] - xp[j])
    return slope * (x - xp[j]) + fp[j]
end

"""
    diel_factor(yday, seconds_of_day)

Instantaneous over daily-mean top-of-atmosphere radiation at `LATITUDE`, in solar time.
"""
function diel_factor(yday, seconds_of_day)
    phi = deg2rad(LATITUDE)
    decl = deg2rad(23.44) * sin(2 * pi * (284 + yday) / 365.0)
    a, b = sin(phi) * sin(decl), cos(phi) * cos(decl)
    h0 = acos(max(-1.0, min(1.0, -a / b)))
    daily_mean = (h0 * a + b * sin(h0)) / pi
    daily_mean <= 0.0 && return 0.0
    hour_angle = 2 * pi * (seconds_of_day / 86400.0 - 0.5)
    return max(0.0, a + b * cos(hour_angle)) / daily_mean
end

"""
    environment(forcing, t_forcing, start, t)

Everything the box is forced with at `t` seconds after `start`. Daily values apply at 12:00 and are
interpolated linearly in time.
"""
function environment(forcing, t_forcing, start, t)
    now = start + Millisecond(round(Int, 1000t))
    yday = Dates.value(now - DateTime(year(now))) / 1000 / 86400.0
    seconds_of_day = Dates.value(now - floor(now, Day)) / 1000

    at(column) = numpy_interp(t, t_forcing, getproperty(forcing, column))

    temp, salt = at(:temp), at(:salt)
    swr = at(:swr)
    DIEL_LIGHT && (swr *= diel_factor(yday, seconds_of_day))

    return (; temp, salt,
              density = 1000.0 + 0.78 * salt - 0.0065 * (temp - 4.0)^2,  # simple approximation (kg/m3)
              swr, wind = at(:wind), yday, Q_river = at(:Q_river))
end

"""
    read_boundary(path)

River and outer-fjord concentrations for every interior state variable, in model units; a variable
the file does not list is 0 (pure dilution).
"""
function read_boundary(path)
    config = YAML.load_file(path)

    return map(("river", "ocean")) do side
        given = config[side]
        for name in keys(given)
            Symbol(name) in INTERIOR || error("$(basename(path)): $side/$name is not an interior state variable")
        end
        NamedTuple{INTERIOR}(map(name -> Float64(get(given, String(name), 0.0)) * get(TO_MODEL_UNITS, String(name), 1.0), INTERIOR))
    end
end

#####
##### The model
#####

"""
    box_rates(configuration, state, env, carbonate)

What `pyfabm.Model.getRates` returns for the box: every state variable's rate of change, the
surface and bottom fluxes divided by the box depth, and the diagnostics.

`carbonate` holds the previous carbonate diagnostics, which the carbonate module falls back on
when its iteration does not converge.
"""
function box_rates(configuration, state, env, carbonate)
    p = configuration.parameters
    light_parameters = configuration.light

    # niva_ecosmo_daylength
    D = day_length(LATITUDE, env.yday)

    # gotm/light at the centre of the single layer
    ext = attenuation_coefficient(p, state.flachl, state.diachl)
    light = gotm_light(env.swr, light_parameters.A, light_parameters.g1, light_parameters.g2,
                       0.0 + DEPTH * 0.5, 0.0 + ext * DEPTH * 0.5)

    # niva_ecosmo_carbonate do
    pressure = 0.5 * DEPTH  # dbar
    carbonate = carbonate_diagnostics(state.dic, state.alk, env.temp, env.salt, env.density, pressure,
                                      carbonate, configuration.pH_scale)

    interior = interior_rates(p, state, env.temp, D, light.par, carbonate.Om_cal)

    surface = (oxy = oxygen_surface_flux(state.oxy, env.temp, env.salt, env.wind),
               dic = carbon_dioxide_surface_flux(state.dic, state.alk, env.temp, env.salt, env.density, pressure,
                                                 env.wind, PCO2_AIR, carbonate.pCO2, configuration.pH_scale))

    bottom = bottom_fluxes(p, state, state, env.temp, carbonate.Om_cal, BOTTOM_STRESS)

    # FABM's get_sources: interior sources, then the surface fluxes over the cell thickness, then
    # the bottom fluxes over it
    interior_total = map(INTERIOR) do name
        rate = interior.rates[name] + get(surface, name, 0.0) / DEPTH
        rate + get(bottom.fluxes, name, 0.0) / DEPTH
    end

    rates = merge(NamedTuple{INTERIOR}(interior_total), bottom.sediment)

    diagnostics = merge(
        NamedTuple{Symbol.("CO2_" .* string.(keys(carbonate)))}(values(carbonate)),
        NamedTuple{Symbol.("ECO_" .* string.(keys(interior.diagnostics)))}(values(interior.diagnostics)),
        (light_swr = light.swr,
         light_par = light.par,
         light_par0 = env.swr * (1.0 - light_parameters.A),
         attenuation_coefficient_of_photosynthetic_radiative_flux = ext,
         Chl_s = state.flachl + state.diachl,
         daylength_day_length = D,
         CO2_wind = env.wind,
         ECO_botstrss = BOTTOM_STRESS))

    return (; rates, diagnostics, carbonate)
end

"""
    clip(state)

FABM's `check_state(repair = .true.)`: every state variable clipped to its minimum.
"""
clip(state) = merge(map(max, ECOSMO_MINIMA, state[INTERIOR]), map(x -> max(0.0, x), state[SEDIMENT]))

closed_parameters(p) =
    EcosmoParameters{Float64}(map(name -> name in BURIAL_PARAMETERS ? 0.0 : getfield(p, name), ECOSMO_PARAMETER_NAMES)...)

total_phosphorus(state) = DEPTH * sum(name -> state[name], P_INTERIOR) + sum(name -> state[name], P_BOTTOM)

"""
    run_box(; yaml, forcing, boundary, closed = false)

Integrate the box over the forcing period. Returns the daily output: the state at the start of each
day, the daily-mean diagnostics and exchange flows, and the forcing.
"""
function run_box(; yaml = joinpath(DATA, "fabm.yaml"), forcing = joinpath(DATA, "forcing.csv"),
                   boundary = joinpath(DATA, "boundary.yaml"), closed = false)

    configuration = fabm_configuration(yaml)
    if closed
        configuration = merge(configuration, (; parameters = closed_parameters(configuration.parameters)))
    end

    table = read_forcing(forcing)
    start = DateTime(table.date[1])
    t_forcing = [Dates.value(DateTime(date) - start) / 1000 + 43200.0 for date in table.date]
    n_days = length(table.date)

    river, ocean = read_boundary(boundary)

    state = configuration.initial_values[(INTERIOR..., SEDIMENT...)]
    # FABM's missing values for the carbonate diagnostics, which they hold until first computed
    carbonate = (pH = 0.0, pCO2 = 0.0, CarbA = 0.0, BiCarb = 0.0, Carb = 0.0, Om_cal = 4.0, Om_arg = 3.0)

    p_start = total_phosphorus(state)

    steps_per_day = round(Int, 86400.0 / DT)
    out_state = Vector{typeof(state)}(undef, n_days)
    out_diagnostics = Vector{Any}(undef, n_days)
    out_exchange = zeros(n_days, 2)

    for day in 1:n_days
        out_state[day] = state
        daily = nothing

        for step in 0:steps_per_day-1
            t = ((day - 1) * steps_per_day + step) * DT
            env = environment(table, t_forcing, start, t)

            result = box_rates(configuration, state, env, carbonate)
            carbonate = result.carbonate
            rates = result.rates

            if !closed
                q_river = env.Q_river
                q_ocean = isnothing(Q_OCEAN) ? q_river * env.salt / max(S_OCEAN - env.salt, 0.1) : Q_OCEAN
                exchange = map(INTERIOR) do name
                    c = state[name]
                    (q_river * (river[name] - c) + q_ocean * (ocean[name] - c)) / VOLUME
                end
                rates = merge(rates, NamedTuple{INTERIOR}(map(+, rates[INTERIOR], exchange)))
                out_exchange[day, 1] += q_river / steps_per_day
                out_exchange[day, 2] += q_ocean / steps_per_day
            end

            state = clip(map((c, rate) -> c + DT * rate, state, rates))

            mean_part = map(value -> value / steps_per_day, result.diagnostics)
            daily = isnothing(daily) ? mean_part : map(+, daily, mean_part)
        end

        out_diagnostics[day] = daily
    end

    all(isfinite, state) || error("Non-finite values in the model state")

    if closed
        p_end = total_phosphorus(state)
        println("Total P (mgC/m2): start $p_start, end $p_end, relative change $((p_end - p_start) / p_start)")
        println("(clipping of near-zero pools to their minimum can cause small deviations)")
    end

    dates = [start + Day(n) for n in 0:n_days-1]
    noon_seconds = [Dates.value(date - start) / 1000 for date in dates]
    forcing_out = map(column -> [numpy_interp(t, t_forcing, getproperty(table, column)) for t in noon_seconds],
                      (temp = :temp, salt = :salt, swr = :swr, wind = :wind))

    return (; dates, state = out_state, diagnostics = out_diagnostics, exchange = out_exchange,
              forcing = forcing_out, closed, yaml)
end

#####
##### Output
#####

function write_output(path, output)
    rm(path; force = true)
    NCDataset(path, "c") do ds
        defDim(ds, "time", length(output.dates))
        defVar(ds, "time", output.dates, ("time",))

        for name in keys(output.state[1])
            defVar(ds, output_name(name), [s[name] for s in output.state], ("time",))
        end
        for name in keys(output.diagnostics[1])
            defVar(ds, String(name), [d[name] for d in output.diagnostics], ("time",))
        end

        defVar(ds, "Q_river", output.exchange[:, 1], ("time",); attrib = ["units" => "m3/s"])
        defVar(ds, "Q_ocean", output.exchange[:, 2], ("time",); attrib = ["units" => "m3/s"])
        for (name, values) in pairs(output.forcing)
            defVar(ds, "forcing_$name", values, ("time",))
        end

        ds.attrib["depth_m"] = DEPTH
        ds.attrib["area_m2"] = AREA
        ds.attrib["latitude"] = LATITUDE
        ds.attrib["closed"] = Int(output.closed)
        ds.attrib["fabm_yaml"] = output.yaml
    end
    println("Wrote $path")
    return path
end

function plot_output(path, output)
    days = Dates.value.(Date.(output.dates) .- Date(output.dates[1]))
    state(name) = [s[name] for s in output.state]
    diagnostic(name) = [d[name] for d in output.diagnostics]

    fig = M.Figure(size = (1200, 1000))
    panels = [
        ("Forcing", [("T (degC)", output.forcing.temp), ("S (PSU)", output.forcing.salt), ("SWR/10 (W/m2)", output.forcing.swr ./ 10)]),
        ("Nutrients (umol/L)", [("NO3", state(:no3) ./ TO_MODEL_UNITS["no3"]), ("NH4", state(:nh4) ./ TO_MODEL_UNITS["nh4"]),
                                ("PO4 x10", state(:pho) ./ TO_MODEL_UNITS["pho"] .* 10), ("Si", state(:sil) ./ TO_MODEL_UNITS["sil"])]),
        ("Chlorophyll a (mg/m3)", [("total", state(:flachl) .+ state(:diachl)), ("small (fla)", state(:flachl)), ("large (dia)", state(:diachl))]),
        ("Plankton and detritus (mgC/m3)", [(String(name), state(name)) for name in (:fla, :dia, :microzoo, :mesozoo, :det)]),
        ("Oxygen", [("O2 (mmol/m3)", state(:oxy))]),
        ("pH", [("pH", diagnostic(:CO2_pH))]),
    ]
    for (n, (title, lines)) in enumerate(panels)
        ax = M.Axis(fig[cld(n, 2), mod1(n, 2)]; title, xlabel = "day of $(year(output.dates[1]))")
        for (label, values) in lines
            M.lines!(ax, days, values; label)
        end
        M.axislegend(ax; labelsize = 10)
    end
    M.save(path, fig)
    println("Wrote $path")
    return path
end

function main(args)
    closed = "--closed" in args
    value(flag, default) = (i = findfirst(==(flag), args); isnothing(i) ? default : args[i+1])

    output = run_box(; yaml = value("--yaml", joinpath(DATA, "fabm.yaml")),
                       forcing = value("--forcing", joinpath(DATA, "forcing.csv")),
                       boundary = value("--boundary", joinpath(DATA, "boundary.yaml")),
                       closed)

    out = value("--out", joinpath(HERE, closed ? "output_closed.nc" : "output.nc"))
    write_output(out, output)
    "--no-plot" in args || plot_output(replace(out, r"\.nc$" => ".png"), output)
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
