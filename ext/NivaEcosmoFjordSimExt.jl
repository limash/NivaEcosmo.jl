"""
ECOSMO in a FjordSim setup: the model, its air-sea gas exchange, the ECOSMO tracers on the open
boundaries, in the initial state and in the rivers, and a surface writer. Loaded when a setup loads
both NivaEcosmo and FjordSim (and NCDatasets, which the prepared files are written with).

The tracer values a setup gives are functions, so one mechanism serves a setup with no
biogeochemical data (seasonal profiles) and one reading a reanalysis (TOPAZ5-ECOSMO): each tracer is
`(z, date, source) -> value`, `z` the depth (m, negative down), `date` a `DateTime` and `source` a
`NamedTuple` of the values at that point of the base file's `T` and `S` and of the source dataset's
prepared variables.
"""
module NivaEcosmoFjordSimExt

using Dates: DateTime
using NCDatasets: NCDataset, defVar, dimnames
using Oceananigans: Center, JLD2Writer, KernelFunctionOperation, SpecifiedTimes, TimeInterval, fields, set!
using Oceananigans.BoundaryConditions: getbc

using NivaEcosmo
using NivaEcosmo: ECOSMO, ECOSMO_TRACERS, CARBONATE_FIELDS, QuadraticBottomStress, ecosmo_sediment,
    ecosmo_surface_fluxes

using FjordSim: BoundaryConditions, BoundaryData, InitialState, Rivers, Simulations, Utils
using FjordSim.Configs:
    AbstractBiogeochemistryConfig,
    AbstractBoundaryConditionConfig,
    AbstractBoundaryDataConfig,
    AbstractInitialStateConfig,
    AbstractRiverConfig,
    AbstractRunConfig,
    AbstractWriterConfig

# The prepared files' own deflate level, so appended variables compress like the ones beside them.
const DEFLATE_LEVEL = 5
# Boundary records computed and written at a time: a month of hourly records.
const BOUNDARY_RECORDS_PER_WRITE = 744
# The base file's own variables a tracer function also sees in `source`, beside the source dataset's.
const BASE_INPUTS = (:T, :S)

# --- The model --------------------------------------------------------------------------------

struct EcosmoBiogeochemistryConfig{C,S} <: AbstractBiogeochemistryConfig
    configuration::C
    surface_shortwave::S
    bottom_drag_coefficient::Float64
    start_day_of_year::Float64
end

NivaEcosmo.ecosmo_biogeochemistry_config(; configuration, surface_shortwave, bottom_drag_coefficient,
                                         start_day_of_year) =
    EcosmoBiogeochemistryConfig(configuration, surface_shortwave, Float64(bottom_drag_coefficient),
                                Float64(start_day_of_year))

# The sediment pools are set here rather than through `initial_conditions`: they are not model fields.
function Simulations.model_biogeochemistry(config::EcosmoBiogeochemistryConfig, grid)
    (; parameters, pH_scale, initial_values) = config.configuration
    bottom_stress = QuadraticBottomStress(drag_coefficient = config.bottom_drag_coefficient)
    sediment = ecosmo_sediment(grid; parameters, bottom_stress)
    for name in (:sed1, :sed2, :sed3, :sed4)
        set!(sediment.fields[name], initial_values[name])
    end
    return ECOSMO(grid; parameters, surface_shortwave = config.surface_shortwave, sediment,
                  start_day_of_year = config.start_day_of_year, pH_scale)
end

# --- Air-sea gas exchange -----------------------------------------------------------------------

struct EcosmoGasExchange <: AbstractBoundaryConditionConfig
    wind_speed::Float64
    atmospheric_pCO2::Float64
    pH_scale::Int
end

NivaEcosmo.ecosmo_gas_exchange(; wind_speed, atmospheric_pCO2, pH_scale) =
    EcosmoGasExchange(Float64(wind_speed), Float64(atmospheric_pCO2), Int(pH_scale))

BoundaryConditions.boundary_condition_sides(config::EcosmoGasExchange, grid, forcing, boundary_data_config, tracers) =
    map(condition -> (top = condition,),
        ecosmo_surface_fluxes(eltype(grid); wind_speed = config.wind_speed,
                              atmospheric_pCO2 = config.atmospheric_pCO2, pH_scale = config.pH_scale))

# --- Open boundaries ----------------------------------------------------------------------------

struct EcosmoBoundaries{B,S,T} <: AbstractBoundaryDataConfig
    base::B
    sources::S
    tracers::T
    data_root::String
    output_file::String
    plot_file::String
    open_edges::Vector{Symbol}
    start_date::DateTime
    stop_date::DateTime
end

function NivaEcosmo.ecosmo_boundaries(base; sources = nothing, tracers)
    isnothing(sources) || sources.open_edges == base.open_edges ||
        throw(ArgumentError("`sources` must open the edges `base` does, $(base.open_edges)."))
    return EcosmoBoundaries(base, sources, tracers, base.data_root, base.output_file, base.plot_file,
                            base.open_edges, base.start_date, base.stop_date)
end

BoundaryData.boundary_variable_names(config::EcosmoBoundaries) =
    merge(BoundaryData.boundary_variable_names(config.base),
          Dict(String(name) => String(name) for name in keys(config.tracers)))

function BoundaryData.fetch_boundary_data(target_grid, config::EcosmoBoundaries)
    downloaded = BoundaryData.fetch_boundary_data(target_grid, config.base)
    BoundaryData.fetch_boundary_data(target_grid, config.sources)
    return downloaded
end

function BoundaryData.regrid_boundary_data(target_grid, config::EcosmoBoundaries)
    result = BoundaryData.regrid_boundary_data(target_grid, config.base)
    sources = BoundaryData.regrid_boundary_data(target_grid, config.sources)
    source_file = isnothing(sources) ? nothing : sources.output_file
    appended = append_boundary_tracers!(result.output_file, config.open_edges, config.tracers, source_file)
    return (; result..., variables = [result.variables; appended])
end

"""
    append_boundary_tracers!(filepath, edges, tracers, source_file)

Write `<edge>_<tracer>` for every tracer and open edge into the prepared boundary file `filepath`,
evaluated at each boundary cell, level and record. `source_file`, a prepared boundary file of the
source dataset (or `nothing`), supplies each point's `source` values, interpolated linearly in time
onto `filepath`'s records, which its axis must span.
"""
function append_boundary_tracers!(filepath, edges, tracers, source_file)
    appended = String[]
    NCDataset(filepath, "a") do ds
        z = Float64.(ds["Nz"][:])
        dates = DateTime.(ds["time"][:])
        with_source(source_file) do source
            records = source_records(source, dates, filepath)
            for edge in edges
                along = edge in (:south, :north) ? "Nx" : "Ny"
                count = ds.dim[along]
                own = filter(pair -> first(pair) in BASE_INPUTS, source_variables(ds, "$(edge)_", (along, "Nz", "time")))
                inputs = source_variables(source, "$(edge)_", (along, "Nz", "time"))
                variables = map(collect(pairs(tracers))) do (name, _)
                    variable_name = BoundaryData.boundary_variable_name(edge, String(name))
                    push!(appended, variable_name)
                    @info "Writing $variable_name into $(basename(filepath))"
                    defVar(ds, variable_name, Float32, (along, "Nz", "time");
                           chunksizes = [count, length(z), 1], deflatelevel = DEFLATE_LEVEL,
                           attrib = ["_FillValue" => NaN32])
                end
                for block in Iterators.partition(eachindex(dates), BOUNDARY_RECORDS_PER_WRITE)
                    block_inputs = merge(read_block(ds, own, [(t, t, 0.0) for t in block]),
                                         read_block(source, inputs, records[block]))
                    for (variable, tracer) in zip(variables, values(tracers))
                        variable[:, :, block] = evaluate(tracer, z, dates[block], block_inputs, count)
                    end
                end
            end
        end
    end
    return appended
end

# `tracer` at each (point, level, record) of a block of boundary records.
function evaluate(tracer, z, dates, inputs, count)
    output = Array{Float32}(undef, count, length(z), length(dates))
    for t in eachindex(dates), k in eachindex(z), i in 1:count
        output[i, k, t] = tracer(z[k], dates[t], map(values -> values[i, k, t], inputs))
    end
    return output
end

with_source(body, ::Nothing) = body(nothing)
with_source(body, filepath) = NCDataset(body, filepath)

source_records(::Nothing, dates, filepath) = [(index, index, 0.0) for index in eachindex(dates)]

# For each of `dates`, the source records bracketing it and the weight of the later one: the source
# (TOPAZ's daily means) is interpolated linearly in time onto the boundary file's hourly axis.
function source_records(source, dates, filepath)
    times = DateTime.(source["time"][:])
    return map(dates) do date
        upper = findfirst(>=(date), times)
        (isnothing(upper) || (upper == 1 && times[1] > date)) && error(
            "The source's prepared boundary file runs $(first(times)) to $(last(times)), which does " *
            "not span $date, a record of $filepath. Prepare it over the same period.",
        )
        times[upper] == date && return (upper, upper, 0.0)
        lower = upper - 1
        (lower, upper, (date - times[lower]) / (times[upper] - times[lower]))
    end
end

"""
    source_variables(source, prefix, dimensions)

The names, without `prefix`, of the source file's variables on `dimensions` whose names start with
`prefix`, as `Symbol => file name` pairs.
"""
source_variables(::Nothing, prefix, dimensions) = Pair{Symbol,String}[]
source_variables(source, prefix, dimensions) = [
    Symbol(chopprefix(name, prefix)) => name
    for name in keys(source) if startswith(name, prefix) && dimnames(source[name]) == dimensions
]

read_block(::Nothing, inputs, records) = (;)
function read_block(source, inputs, records)
    span = minimum(first, records):maximum(record -> record[2], records)
    offset = first(span) - 1
    return NamedTuple(map(inputs) do (name, file)
        data = Float32.(coalesce.(source[file][:, :, span], NaN32))
        blended = similar(data, size(data, 1), size(data, 2), length(records))
        for (t, (lower, upper, weight)) in enumerate(records)
            @views @. blended[:, :, t] = (1 - weight) * data[:, :, lower - offset] + weight * data[:, :, upper - offset]
        end
        name => blended
    end)
end

# --- Initial state ------------------------------------------------------------------------------

struct EcosmoInitialState{B,S,T} <: AbstractInitialStateConfig
    base::B
    sources::S
    tracers::T
    data_root::String
    output_file::String
    start_date::DateTime
end

NivaEcosmo.ecosmo_initial_state(base; sources = nothing, tracers) =
    EcosmoInitialState(base, sources, tracers, base.data_root, base.output_file, base.start_date)

function InitialState.fetch_initial_state(target_grid, config::EcosmoInitialState, date)
    downloaded = InitialState.fetch_initial_state(target_grid, config.base, date)
    isnothing(config.sources) || InitialState.fetch_initial_state(target_grid, config.sources, date)
    return downloaded
end

function InitialState.regrid_initial_state(target_grid, config::EcosmoInitialState, date; edges = nothing)
    result = InitialState.regrid_initial_state(target_grid, config.base, date; edges)
    source_file = isnothing(config.sources) ? nothing :
                  InitialState.regrid_initial_state(target_grid, config.sources, date; edges).output_file
    append_initial_tracers!(result.output_file, config.tracers, source_file, date)
    return (; result..., variables = [result.variables; map(String, collect(keys(config.tracers)))])
end

"""
    append_initial_tracers!(filepath, tracers, source_file, date)

Write every tracer into the prepared initial state `filepath`, evaluated at each water cell (where
its `T` is finite) at `date`, with `source` from the source dataset's prepared initial state.
"""
function append_initial_tracers!(filepath, tracers, source_file, date)
    NCDataset(filepath, "a") do ds
        z = Float64.(ds["Nz"][:])
        water = isfinite.(coalesce.(ds["T"][:, :, :, 1], NaN32))
        dimensions = ("Nx", "Ny", "Nz", "time")
        state(file, variables) = NamedTuple(name => Float32.(coalesce.(file[filename][:, :, :, 1], NaN32))
                                            for (name, filename) in variables)
        own = filter(pair -> first(pair) in BASE_INPUTS, source_variables(ds, "", dimensions))
        inputs = merge(state(ds, own),
                       with_source(source -> state(source, source_variables(source, "", dimensions)), source_file))
        for (name, tracer) in pairs(tracers)
            @info "Writing $name into $(basename(filepath))"
            values = fill(NaN32, size(water))
            for index in CartesianIndices(water)
                water[index] || continue
                values[index] = tracer(z[index[3]], date, map(source -> source[index], inputs))
            end
            variable = defVar(ds, String(name), Float32, ("Nx", "Ny", "Nz", "time");
                              deflatelevel = DEFLATE_LEVEL, attrib = ["_FillValue" => NaN32])
            variable[:, :, :, 1] = values
        end
    end
    return filepath
end

# --- Rivers -------------------------------------------------------------------------------------

struct EcosmoRivers{R,C} <: AbstractRiverConfig
    base::R
    concentrations::C
    data_root::String
    output_file::String
    plot_file::String
    search_radius::Int
    maximum_inflow_velocity::Float64
end

NivaEcosmo.ecosmo_rivers(base, concentrations) =
    EcosmoRivers(base, concentrations, base.data_root, base.output_file, base.plot_file,
                 base.search_radius, base.maximum_inflow_velocity)

Utils.coverage_window(rivers::EcosmoRivers) = Utils.coverage_window(rivers.base)
Rivers.river_locations(rivers::EcosmoRivers) = Rivers.river_locations(rivers.base)
Rivers.river_discharge(rivers::EcosmoRivers, times) = Rivers.river_discharge(rivers.base, times)
Rivers.fetch_rivers(target_grid, rivers::EcosmoRivers) = Rivers.fetch_rivers(target_grid, rivers.base)

function Rivers.river_series(rivers::EcosmoRivers, times)
    series = Rivers.river_series(rivers.base, times)
    count = length(Rivers.river_locations(rivers))
    for (name, concentration) in pairs(rivers.concentrations)
        series[String(name)] = Float32[concentration(time) for _ in 1:count, time in times]
    end
    return series
end

# --- Surface writer -----------------------------------------------------------------------------

struct EcosmoSurfaceWriter <: AbstractWriterConfig
    name::Symbol
    output_file::String
    interval::Float64
    split_times::Vector{Float64}
    overwrite_existing::Bool
end

NivaEcosmo.ecosmo_surface_writer(; name, output_file, interval, split_times, overwrite_existing) =
    EcosmoSurfaceWriter(Symbol(name), String(output_file), Float64(interval), Float64.(split_times),
                        Bool(overwrite_existing))

Simulations.output_path_trait(::EcosmoSurfaceWriter) = Simulations.NamesOutputFile()

# Oceananigans' fluxes are positive upward; FABM's, written here, into the water.
@inline air_sea_flux(i, j, k, grid, condition, clock, model_fields) =
    -getbc(condition, i, j, grid, clock, model_fields)

function surface_outputs(model)
    model_fields = fields(model)
    Nz = size(model.grid, 3)
    names = (:T, :S, :u, :v, :e, ECOSMO_TRACERS..., CARBONATE_FIELDS..., :PAR)
    surface = NamedTuple{names}(map(name -> view(model_fields[name], :, :, Nz), names))
    flux(name) = KernelFunctionOperation{Center,Center,Nothing}(
        air_sea_flux, model.grid, model.tracers[name].boundary_conditions.top, model.clock, model_fields)
    return merge(surface, (η = model_fields.η, day_length = model_fields.day_length),
                 model.biogeochemistry.sediment.fields,
                 (oxy_surface_flux = flux(:oxy), dic_surface_flux = flux(:dic)))
end

function Simulations.attach_writer!(simulation, writer::EcosmoSurfaceWriter, config::AbstractRunConfig, loop)
    ocean = simulation.model.ocean
    haskey(ocean.output_writers, writer.name) && close(pop!(ocean.output_writers, writer.name))
    filepath = Simulations.loop_output_path(writer, config, loop)
    mkpath(dirname(filepath))
    ocean.output_writers[writer.name] = JLD2Writer(
        ocean.model, surface_outputs(ocean.model);
        filename = filepath, schedule = TimeInterval(writer.interval),
        file_splitting = SpecifiedTimes(writer.split_times),
        overwrite_files = writer.overwrite_existing, with_halos = false,
    )
    return simulation
end

end # module
