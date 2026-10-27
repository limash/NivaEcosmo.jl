# Møre og Romsdal: the coast from Stad to Hitra (5.12-9.00°E, 62.00-63.43°N), with Storfjorden,
# Romsdalsfjorden, the Kristiansund fjords, the Møre shelf and the continental slope off it in the
# north-western corner. A NorKyst-800m initial state, hourly NorKyst-800m data on the open western
# and northern edges, NVE rivers, NORA3 atmosphere and the NIVA ECOSMO biogeochemistry
# (NivaEcosmo.jl), run over 2023.
#
# # Preparing and running it
#
# Installation and the order of the steps are in `fjordsim_model/README.md`. Every step is
#
#   julia --project=fjordsim_model -m FjordSim <step> --config fjordsim_model/moro.jl
#
# `add_rivers` needs `NVE_API_KEY` (free at https://hydapi.nve.no/Users). Data goes to
# `~/FjordSim_data/moro/`, the Geonorge FileGDB to its own directory there, `geonorge/`, downloaded
# (4.5 GB) by `prepare_bathymetry` if it is missing. Results go to `~/FjordSim_results/moro/`.
#
# # What follows `examples/oslofjorden_validation.jl`, and what does not
#
# The physics is the validation run's, with every coefficient that depends on grid spacing scaled
# from 193 m to 500 m: bounded WENO for salinity, CATKE's `minimum_tke` at 7e-7, the biharmonic
# tracer diffusivity κ at 0, `oslofjorden.jl`'s levels with 5.06 m levels between 41 and 132 m,
# its bathymetry knobs, sponge, bottom drag, open-boundary timescales and time stepping. Each is
# argued beside its value. Three things do not carry over. The validation run reads the Norkyst-v3
# hindcast because its window predates 2017, while 2023 is in the operational archive. Its tide
# correction was fitted at Viker. Its initial state is a previous run's.
#
# # The biogeochemistry
#
# ECOSMO's 17 tracers ride beside T and S, with the A20 v3aa parameters of NivaEcosmo's
# `data/fabm.yaml`. As in `examples/oslofjorden_npzd.jl`, FjordSim's open boundaries need a
# prepared series for every tracer and the rivers may carry any of them, so this file adds three
# things:
#
# - `EcosmoBoundaries` appends the tracers' profiles to the prepared boundary file;
# - `WithProfiles` starts them from the same profiles, beside NorKyst's T, S, u and v;
# - `EcosmoRivers` gives every river one concentration per tracer.
#
# **Every biogeochemical value is a placeholder, not a measurement.** They are in
# `data/moro_ecosmo.yaml`: typical values for Atlantic-influenced Norwegian coastal water and
# west-coast rivers, with a seasonal cycle. On the boundaries the surface nutrients are drawn down
# and the plankton grows through a spring-to-autumn window, uniform along the edges. In the rivers
# each concentration is a cosine in day of year. Replace them with observations or a reanalysis
# before the run is used for anything but bringing the setup up. Light is a clear-sky shortwave
# scaled by a constant clearness index, and the gas exchange runs at a constant wind. NORA3 has
# both, but reading them here would mean reaching into NumericalEarth's atmosphere objects.

using FjordSim
using FjordSim: BoundaryConditions, BoundaryData, Rivers, Simulations
using FjordSim.Bathymetry: GEONORGE_DYBDEDATA_GDB
using FjordSim.Configs:
    AbstractBiogeochemistryConfig,
    AbstractBoundaryConditionConfig,
    AbstractBoundaryDataConfig,
    AbstractRiverConfig,
    AbstractWriterConfig,
    run_config,
    start_date
using Oceananigans:
    CATKEVerticalDiffusivity,
    Center,
    HydrostaticSphericalCoriolis,
    JLD2Writer,
    KernelFunctionOperation,
    SeawaterBuoyancy,
    SpecifiedTimes,
    TimeInterval,
    WENO,
    WENOVectorInvariant,
    fields,
    set!
using Oceananigans.BoundaryConditions: getbc
using Oceananigans.TurbulenceClosures: HorizontalScalarBiharmonicDiffusivity
using Oceananigans.Units: day, days, hour, hours, minutes, second
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using NumericalEarth: FreezingLimitedOceanTemperature
using NivaEcosmo:
    NivaEcosmo,
    ECOSMO,
    ECOSMO_TRACERS,
    QuadraticBottomStress,
    ecosmo_sediment,
    ecosmo_surface_fluxes,
    fabm_configuration
using NCDatasets: NCDataset, defVar
using YAML
using Dates: DateTime, Millisecond, Month, dayofyear

data_root = fjord_data_root("moro")
FT = Float32
years = [2023]

# Stated once because two things read it: the run's clock, and the day of the year that ECOSMO's
# day length and the clear-sky shortwave count from.
const START_DATE = DateTime(2023, 1, 1)
# Days since midnight on 1 January, as ECOSMO's `start_day_of_year` counts them.
const START_DAY_OF_YEAR = dayofyear(START_DATE) - 1
# Midnight on the first of each later month of the year, in seconds of model time: where the hourly
# surface writer starts a new file.
const MONTH_STARTS = [(START_DATE + Month(month) - START_DATE) / Millisecond(1000) for month in 1:11]

# Shared by the momentum drag and the sediment's resuspension stress, as NivaEcosmo asks: the
# sediment's bottom stress should come from the same drag law as the ocean's.
const BOTTOM_DRAG_COEFFICIENT = 0.003

# --- The biogeochemical placeholders -----------------------------------------------------------

# Every number is in `data/moro_ecosmo.yaml`, which says why each has the value it has.
const MORO_ECOSMO = YAML.load_file(joinpath(pkgdir(NivaEcosmo), "data", "moro_ecosmo.yaml"))

# The file's units to ECOSMO's. N, P and Si become mgC m⁻³ through the Redfield ratios the model
# itself uses.
const TO_ECOSMO_UNITS = Dict(
    "umol N/L" => NivaEcosmo.REDFIELD_C_N * NivaEcosmo.CARBON_MOLAR_MASS,
    "umol P/L" => NivaEcosmo.REDFIELD_C_P * NivaEcosmo.CARBON_MOLAR_MASS,
    "umol Si/L" => NivaEcosmo.REDFIELD_C_SI * NivaEcosmo.CARBON_MOLAR_MASS,
    "mmol/m3" => 1.0,
    "mgC/m3" => 1.0,
)

const SEASON = MORO_ECOSMO["season"]
const BOUNDARY = MORO_ECOSMO["boundary"]

"""
    seasonal_window(date)

How far into the productive season `date` is, from 0 in midwinter to 1 through the summer: a
logistic opening at `onset_day` times one closing at `end_day`, as `oslofjorden_npzd.jl`'s
`seasonal_drawdown`.
"""
function seasonal_window(date)
    day = dayofyear(date)
    onset = 1 / (1 + exp(-(day - SEASON["onset_day"]) / SEASON["onset_width"]))
    ending = 1 / (1 + exp((day - SEASON["end_day"]) / SEASON["end_width"]))
    return onset * ending
end

"""
    water_mass_profile(z, surface, deep)

`surface` above `surface_layer_depth`, `deep` below `deep_water_depth`, linear between, at depth
`z` (negative down, metres).
"""
function water_mass_profile(z, surface, deep)
    top = BOUNDARY["surface_layer_depth"]
    bottom = BOUNDARY["deep_water_depth"]
    weight = clamp((-z - top) / (bottom - top), 0, 1)
    return (1 - weight) * surface + weight * deep
end

"""
    dissolved_profile(entry)

A dissolved tracer's `(z, date) -> value`, in ECOSMO's units: the winter column between the
surface and deep values, plus the summer change at the surface times `seasonal_window`, decaying
with depth over `depletion_scale`.
"""
function dissolved_profile(entry)
    factor = TO_ECOSMO_UNITS[entry["unit"]]
    winter, summer, deep = entry["winter"], entry["summer"], entry["deep"]
    scale = BOUNDARY["depletion_scale"]
    return (z, date) -> factor * (water_mass_profile(z, winter, deep) +
                                  (summer - winter) * seasonal_window(date) * exp(z / scale))
end

"""
    biological_profile(entry)

Plankton's or organic matter's `(z, date) -> value`, in ECOSMO's units: a surface value moving
from winter to summer with `seasonal_window`, decaying with depth over `euphotic_scale` to the deep
value.
"""
function biological_profile(entry)
    factor = TO_ECOSMO_UNITS[entry["unit"]]
    winter, summer, deep = entry["winter"], entry["summer"], entry["deep"]
    scale = BOUNDARY["euphotic_scale"]
    return (z, date) -> factor * (deep + (winter + (summer - winter) * seasonal_window(date) - deep) * exp(z / scale))
end

"""
    ecosmo_profiles()

Every ECOSMO tracer as `(z, date) -> value`, in `ECOSMO_TRACERS` order. A tracer the file does not
give is a `KeyError` here, rather than a missing boundary series at the start of a run.
"""
function ecosmo_profiles()
    profiles = Dict{Symbol,Function}()
    for (name, entry) in BOUNDARY["dissolved"]
        profiles[Symbol(name)] = dissolved_profile(entry)
    end
    for (name, entry) in BOUNDARY["biological"]
        profiles[Symbol(name)] = biological_profile(entry)
    end
    for (name, entry) in BOUNDARY["chlorophyll"]
        carbon = profiles[Symbol(entry["from"])]
        ratio = entry["carbon_to_chlorophyll"]
        profiles[Symbol(name)] = (z, date) -> carbon(z, date) / ratio
    end
    return NamedTuple{ECOSMO_TRACERS}(map(name -> profiles[name], ECOSMO_TRACERS))
end

# One `NamedTuple` because the initial conditions and the boundary append need the same set,
# evaluated the same way.
const ECOSMO_PROFILES = ecosmo_profiles()

"""
    river_concentrations()

Each tracer the file's `river` section names, as `date -> value` in ECOSMO's units: a cosine in
day of year around its mean, peaking at `peak_day`.
"""
function river_concentrations()
    concentrations = map(collect(MORO_ECOSMO["river"])) do (name, entry)
        factor = TO_ECOSMO_UNITS[entry["unit"]]
        mean, amplitude, peak = entry["mean"], entry["amplitude"], entry["peak_day"]
        Symbol(name) => date -> factor * (mean + amplitude * cos(2π * (dayofyear(date) - peak) / 365))
    end
    return (; concentrations...)
end

# --- Light and wind ----------------------------------------------------------------------------

# `oslofjorden_npzd.jl`'s clear-sky shortwave, as shortwave rather than PAR: `GOTMLight` takes the
# shortwave and keeps its own visible fraction, 1 - A = 0.42. The clearness index is that file's
# value for Oslo, so the cloudier west coast likely gets too much light.
const SOLAR_CONSTANT = 1361.0
const CLEARNESS_INDEX = 0.5

# 10 m wind for the O₂ and CO₂ exchange, m s⁻¹. Assumed. The exchange goes as the square of the
# wind, so a mean speed underestimates it; NORA3's root-mean-square speed over the box would be
# the better constant.
const WIND_SPEED = 7.0

# CO₂ in air, ppm: the 2023 global mean to within about 1 ppm.
const ATMOSPHERIC_PCO2 = 420.0

"""
    surface_shortwave(λ, φ, t)

Daily-mean downwelling shortwave just below the surface, W m⁻², at longitude `λ` and latitude `φ`
(degrees) and `t` seconds after `START_DATE`.

`oslofjorden_npzd.jl`'s `surface_par` without its PAR fraction: Cooper's (1969) declination and the
daily-mean top-of-atmosphere integral over the hour angle, scaled by `CLEARNESS_INDEX`. Evaluated
inside a kernel, so no branches and no allocation.
"""
@inline function surface_shortwave(λ, φ, t)
    day_of_year = START_DAY_OF_YEAR + t / 86400
    declination = deg2rad(23.44) * sin(2π * (day_of_year + 284) / 365)
    latitude = deg2rad(φ)

    hour_angle = acos(clamp(-tan(latitude) * tan(declination), -1, 1))
    eccentricity = 1 + 0.033 * cos(2π * day_of_year / 365)

    insolation =
        SOLAR_CONSTANT / π * eccentricity * (
            hour_angle * sin(latitude) * sin(declination) +
            cos(latitude) * cos(declination) * sin(hour_angle)
        )

    return CLEARNESS_INDEX * max(0, insolation)
end

# --- ECOSMO as a FjordSim biogeochemistry ------------------------------------------------------

"""
    EcosmoConfig(; configuration, surface_shortwave, bottom_drag_coefficient, start_day_of_year)

The knobs `model_biogeochemistry` builds `ECOSMO` from, once the grid exists. `configuration` is
what `fabm_configuration` returns.

The light is `GOTMLight` with `gotm/light`'s defaults (`g2 = 23 m`), not the `fabm.yaml`'s
`g2 = 4 m`, which that file sets for the turbid Drammensfjord and marks as a placeholder.
"""
struct EcosmoConfig{C,S} <: AbstractBiogeochemistryConfig
    configuration::C
    surface_shortwave::S
    bottom_drag_coefficient::Float64
    start_day_of_year::Float64
end

EcosmoConfig(; configuration, surface_shortwave, bottom_drag_coefficient, start_day_of_year) =
    EcosmoConfig(configuration, surface_shortwave, Float64(bottom_drag_coefficient), Float64(start_day_of_year))

"""
    model_biogeochemistry(config::EcosmoConfig, grid)

`ECOSMO` on `grid`, with its sediment pools set to the `fabm.yaml`'s initial values. Set here
rather than through `initial_conditions`, because those are model fields and the sediment's are
not.
"""
function Simulations.model_biogeochemistry(config::EcosmoConfig, grid)
    (; parameters, pH_scale, initial_values) = config.configuration
    bottom_stress = QuadraticBottomStress(drag_coefficient = config.bottom_drag_coefficient)
    sediment = ecosmo_sediment(grid; parameters, bottom_stress)

    for name in (:sed1, :sed2, :sed3, :sed4)
        set!(sediment.fields[name], initial_values[name])
    end

    return ECOSMO(
        grid;
        parameters,
        surface_shortwave = config.surface_shortwave,
        sediment,
        start_day_of_year = config.start_day_of_year,
        pH_scale,
    )
end

"""
    EcosmoGasExchange(; wind_speed, atmospheric_pCO2, pH_scale)

ECOSMO's air-sea O₂ and CO₂ exchange as top conditions on `oxy` and `dic`. A piece of its own,
because `AirSeaFluxes` gives a top condition to T and S only. `pH_scale` must be the one `ECOSMO`
gets.
"""
struct EcosmoGasExchange <: AbstractBoundaryConditionConfig
    wind_speed::Float64
    atmospheric_pCO2::Float64
    pH_scale::Int
end

EcosmoGasExchange(; wind_speed, atmospheric_pCO2, pH_scale) =
    EcosmoGasExchange(Float64(wind_speed), Float64(atmospheric_pCO2), Int(pH_scale))

BoundaryConditions.boundary_condition_sides(config::EcosmoGasExchange, grid, forcing, boundary_data_config, tracers) =
    map(
        condition -> (top = condition,),
        ecosmo_surface_fluxes(
            eltype(grid);
            wind_speed = config.wind_speed,
            atmospheric_pCO2 = config.atmospheric_pCO2,
            pH_scale = config.pH_scale,
        ),
    )

"""
    WithProfiles(base, profiles)

`base`'s initial conditions, plus each entry of `profiles`, a `(z, date) -> value`, at the run's
`start_date`.

`FromInitialState` reads only what the NorKyst file carries (T, S, u and v) and leaves every other
tracer at zero, which ECOSMO's state repair would then clip to its minimum.
"""
struct WithProfiles{B,P}
    base::B
    profiles::P
end

function Simulations.resolve_initial_conditions(initial_conditions::WithProfiles, grid, config, tracers)
    state = Simulations.resolve_initial_conditions(initial_conditions.base, grid, config, tracers)
    date = start_date(run_config(config))
    profiles = map(profile -> (λ, φ, z) -> profile(z, date), initial_conditions.profiles)
    return merge(state, profiles)
end

# --- EcosmoBoundaries: the NorKyst boundary file with the ECOSMO tracers appended ---------------

# FjordSim's own writer settings, so the appended variables compress like the ones beside them.
const BOUNDARY_DEFLATE_LEVEL = 5

"""
    append_boundary_variables!(filepath, edges, profiles)

Add one `<edge>_<name>` variable per entry of `profiles` and per edge to the prepared boundary file
at `filepath`, on its along-edge, depth and time axes. Each profile is `(z, date) -> value`.
Uniform along the edge, and finite under land too, which the boundary reader fills from its
neighbours anyway.
"""
function append_boundary_variables!(filepath, edges, profiles)
    NCDataset(filepath, "a") do ds
        z = ds["Nz"][:]
        dates = ds["time"][:]

        for edge in edges
            along = edge in (:south, :north) ? "Nx" : "Ny"
            n_along = ds.dim[along]

            for (name, profile) in pairs(profiles)
                variable_name = BoundaryData.boundary_variable_name(edge, String(name))
                @info "Writing $variable_name into $(basename(filepath))"

                variable = defVar(
                    ds, variable_name, Float32, (along, "Nz", "time");
                    deflatelevel = BOUNDARY_DEFLATE_LEVEL,
                    attrib = ["_FillValue" => NaN32],
                )
                section = Float32[profile(depth, date) for depth in z, date in dates]
                variable[:, :, :] = repeat(reshape(section, 1, length(z), length(dates)), n_along, 1, 1)
            end
        end
    end

    return filepath
end

"""
    EcosmoBoundaries(base)

`base`'s open-boundary data with the ECOSMO tracers' profiles appended to the file it prepares.
Not optional: `OpenLateralBoundaryFromData` radiates every model tracer towards its prepared series
and raises for a tracer with none.
"""
struct EcosmoBoundaries{B} <: AbstractBoundaryDataConfig
    base::B
    data_root::String
    output_file::String
    plot_file::String
    open_edges::Vector{Symbol}
end

EcosmoBoundaries(base) = EcosmoBoundaries(base, base.data_root, base.output_file, base.plot_file, base.open_edges)

BoundaryData.boundary_variable_names(config::EcosmoBoundaries) = merge(
    BoundaryData.boundary_variable_names(config.base),
    Dict(String(name) => String(name) for name in keys(ECOSMO_PROFILES)),
)

BoundaryData.fetch_boundary_data(target_grid, config::EcosmoBoundaries) =
    BoundaryData.fetch_boundary_data(target_grid, config.base)

function BoundaryData.regrid_boundary_data(target_grid, config::EcosmoBoundaries; coverage = nothing)
    result = BoundaryData.regrid_boundary_data(target_grid, config.base; coverage)
    append_boundary_variables!(result.output_file, config.open_edges, ECOSMO_PROFILES)

    appended = [
        BoundaryData.boundary_variable_name(edge, String(name))
        for edge in config.open_edges for name in keys(ECOSMO_PROFILES)
    ]

    return (; result..., variables = [result.variables; appended])
end

# --- EcosmoRivers: NVE's rivers carrying the ECOSMO tracers -------------------------------------

"""
    EcosmoRivers(base, concentrations)

`base`'s rivers, carrying each tracer in `concentrations` (a `NamedTuple` of `date -> value`) as
well as `base`'s own temperature and salinity. Every hook but `river_series` delegates, so the
mouths, columns and discharge are `base`'s. As `oslofjorden_npzd.jl`'s `NPZDRivers`, but seasonal
where that one is constant: `NVERiversConfig`'s `constants` cannot vary in time.
"""
struct EcosmoRivers{R,C} <: AbstractRiverConfig
    base::R
    concentrations::C
    data_root::String
    output_file::String
    plot_file::String
    search_radius::Int
    maximum_inflow_velocity::Float64
end

EcosmoRivers(base, concentrations) = EcosmoRivers(
    base,
    concentrations,
    base.data_root,
    base.output_file,
    base.plot_file,
    base.search_radius,
    base.maximum_inflow_velocity,
)

Rivers.river_locations(rivers::EcosmoRivers) = Rivers.river_locations(rivers.base)
Rivers.river_discharge(rivers::EcosmoRivers, times) = Rivers.river_discharge(rivers.base, times)
Rivers.river_search_radius(rivers::EcosmoRivers) = Rivers.river_search_radius(rivers.base)
Rivers.river_maximum_inflow_velocity(rivers::EcosmoRivers) = Rivers.river_maximum_inflow_velocity(rivers.base)
Rivers.fetch_rivers(target_grid, rivers::EcosmoRivers) = Rivers.fetch_rivers(target_grid, rivers.base)

function Rivers.river_series(rivers::EcosmoRivers, times)
    series = Rivers.river_series(rivers.base, times)
    n_rivers = length(Rivers.river_locations(rivers))
    for (name, concentration) in pairs(rivers.concentrations)
        series[String(name)] = Float32[concentration(time) for _ in 1:n_rivers, time in times]
    end
    return series
end

# --- EcosmoSurfaceWriter: the surface of everything ECOSMO carries -----------------------------

"""
    EcosmoSurfaceWriter(; name, output_file, interval, split_times, overwrite_existing)

A JLD2 snapshot every `interval` of the top cell of T, S, u, v, e, every ECOSMO tracer, the carbonate
system (`pH`, `pCO2`, `Carb`, `Om_cal`) and `PAR`; of `η` and `day_length`; of the sediment pools
`sed1`-`sed4`, one per column; and of the air-sea fluxes `oxy_surface_flux` (mmol O₂ m⁻² s⁻¹) and
`dic_surface_flux` (mmol C m⁻² s⁻¹), positive into the water as in FABM. A new file part,
`_part<n>.jld2`, starts at each of `split_times`, in seconds of model time.

Its own writer because `FieldSnapshotWriter` writes whole fields, and only those `Oceananigans.fields`
names: the sediment pools are not among them, and the gas fluxes exist only as boundary conditions.
"""
struct EcosmoSurfaceWriter <: AbstractWriterConfig
    name::Symbol
    output_file::String
    interval::Float64
    split_times::Vector{Float64}
    overwrite_existing::Bool
end

EcosmoSurfaceWriter(; name, output_file, interval, split_times, overwrite_existing) = EcosmoSurfaceWriter(
    Symbol(name),
    String(output_file),
    Float64(interval),
    Float64.(split_times),
    Bool(overwrite_existing),
)

Simulations.output_path_trait(::EcosmoSurfaceWriter) = Simulations.NamesOutputFile()

# The boundary conditions are positive upward, as Oceananigans' fluxes are.
@inline air_sea_flux(i, j, k, grid, condition, clock, model_fields) =
    -getbc(condition, i, j, grid, clock, model_fields)

"""
    ecosmo_surface_outputs(model)

What `EcosmoSurfaceWriter` writes, as the `NamedTuple` `JLD2Writer` takes. The 3D fields are views of
their top level, which share the model's memory; the fluxes are computed when a record is written,
with the fields the boundary conditions read.
"""
function ecosmo_surface_outputs(model)
    model_fields = fields(model)
    Nz = size(model.grid, 3)

    names = (:T, :S, :u, :v, :e, ECOSMO_TRACERS..., NivaEcosmo.CARBONATE_FIELDS..., :PAR)
    surface = NamedTuple{names}(map(name -> view(model_fields[name], :, :, Nz), names))

    flux(name) = KernelFunctionOperation{Center,Center,Nothing}(
        air_sea_flux,
        model.grid,
        model.tracers[name].boundary_conditions.top,
        model.clock,
        model_fields,
    )

    return merge(
        surface,
        (η = model_fields.η, day_length = model_fields.day_length),
        model.biogeochemistry.sediment.fields,
        (oxy_surface_flux = flux(:oxy), dic_surface_flux = flux(:dic)),
    )
end

"""
    attach_writer!(simulation, writer::EcosmoSurfaceWriter, config, loop)

Attach a `JLD2Writer` to the ocean simulation, as FjordSim's `FieldSnapshotWriter` does: a fresh one
each loop, with a fresh splitting schedule, since the clock restarts at zero.
"""
function Simulations.attach_writer!(simulation, writer::EcosmoSurfaceWriter, config, loop)
    ocean_simulation = simulation.model.ocean

    haskey(ocean_simulation.output_writers, writer.name) &&
        close(pop!(ocean_simulation.output_writers, writer.name))

    filepath = Simulations.loop_output_path(writer, config, loop)
    mkpath(dirname(filepath))

    ocean_simulation.output_writers[writer.name] = JLD2Writer(
        ocean_simulation.model,
        ecosmo_surface_outputs(ocean_simulation.model);
        filename = filepath,
        schedule = TimeInterval(writer.interval),
        file_splitting = SpecifiedTimes(writer.split_times),
        overwrite_files = writer.overwrite_existing,
        with_halos = false,
    )

    return simulation
end

# --- The config ---------------------------------------------------------------------------------

# The A20 v3aa parameters, unchanged, with the carbonate module's pH scale and the sediment's
# initial pools. Only the light and the initial water column differ from the file, both above. The
# parameters are in the grid's `FT`, which `ECOSMO` requires.
ecosmo_configuration = fabm_configuration(joinpath(pkgdir(NivaEcosmo), "data", "fabm.yaml"), FT)

# `oslofjorden.jl`'s levels with 41-132 m split into eighteen 5.06 m levels, as in
# `examples/oslofjorden_validation.jl`, where the split cut the Drammensfjord's salt-budget
# diffusivity below 41 m by 2.6-7.7 times in 60 days. The mechanism, WENO's implicit vertical
# diffusion ~|w|Δz on thick levels, is not specific to that grid; here it acts on the sill basins of
# the fjords and on the coastal current's lower interface.
#
# Below 132 m, `oslofjorden.jl`'s 33.5 m layers, continued down the continental slope to -1137 m,
# against a deepest cell of 1153 m in `bathymetry.nc`. The 30 columns below the face, by at most
# 16 m, are on the slope within 11 cells of the open edges, inside the sponge, and `PartialCellBottom`
# clamps them to it. Thicker deep layers would cost the fjord basins, which reach 688 m, their
# vertical structure. `grid.Lz` of 1137 m gives a barotropic wave speed of 106 m/s, so the split-explicit substep is
# about 3 s.
#
# 59 levels on 396 x 315 columns is 7.4 M cells, the size of the Oslofjord validation grid.
z_faces = vcat(
    collect(range(-1137.0, -165.5; step = 33.5)),
    collect(range(-132.0, -41.0; length = 19)),
    [-32.0, -25.0, -19.0, -14.5, -10.8, -7.9, -5.5, -3.7, -2.2, -1.0, 0.0],
)

FjordConfig(
    grid_config = EvenGrid(
        # 0.0098° by 0.00454°: 512 m by 505 m at 62.0°N and 488 m by 505 m at 63.43°N, 2.6 times
        # Oslofjord's 193 m. Every coefficient below that depends on grid spacing is scaled from
        # `oslofjorden_validation.jl` by Δx = 500 m, not copied.
        size      = (396, 315, length(z_faces) - 1),
        halo      = (7, 7, 7),
        longitude = (5.12, 9.0),
        latitude  = (62.00, 63.43),
        # `z_faces` is written into `bathymetry.nc` and `simulation_grid` reads the file, so a
        # change here needs `prepare_bathymetry` re-run.
        z_faces   = z_faces,
        float_type = FT,
    ),
    # `oslofjorden.jl`'s knobs, which the validation run takes unchanged.
    bathymetry_config = DybdedataConfig(
        data_root             = data_root,
        output_file           = "bathymetry.nc",
        plot_file             = "bathymetry.png",
        # The national FileGDB, in its own directory under this setup's data root.
        geodatabase_file      = joinpath("geonorge", GEONORGE_DYBDEDATA_GDB),
        raw_resolution_factor = 4,
        padding_cells         = 2,
        include_contours      = false,
        contour_stride        = 10,
        interpolation_passes  = 1,
        major_basins          = 1,
        minimum_depth         = 2.0,
        # Off, unlike Oslofjord's 5: the northern open edge runs through the Smøla and Hitra
        # skerries, and the stage would flood the islands within the edge band.
        open_boundary_land_cells = 0,
        # 6 cells, 1.5 km² at this resolution. Check that the Smøla and Hitra skerries keep their
        # sounds open.
        max_island_cells      = 6,
        close_narrow_passages = true,
        spike_ratio           = 0.5,
        # Must equal the `minimum_fractional_cell_height` the grid gives `PartialCellBottom`.
        minimum_cell_fraction = 0.2,
        max_slope_factor      = 0.5,
        preserve_sills        = true,
        geonorge_cache        = true,
        regrid_cache          = false,
    ),
    # One NorKyst daily mean, the day `start_date` falls in, for `FromInitialState()` below.
    initial_state_config = NorKystInitialStateConfig(
        data_root        = data_root,
        output_directory = "norkyst_initial",
        output_file      = "initial_state.nc",
        parameters       = ["temperature", "salinity", "u_eastward", "v_northward"],
    ),
    # Salinity 0 (`NVERiversConfig`'s default `constants`), NVE's temperature, and the seasonal
    # ECOSMO concentrations of `data/moro_ecosmo.yaml`.
    river_config = EcosmoRivers(
        NVERiversConfig(
            data_root   = data_root,
            output_file = "rivers_nve.nc",
            plot_file   = "rivers_nve.png",
            years       = years,
            # Mouths are discovered from NVE's ELVIS network and REGINE catchments at 0.5 m³/s, as on
            # `oslofjorden.jl`, where that kept 99.6 % of the discharge.
            minimum_discharge = 0.5,
            # The vertical Courant bound on the 1 m top cell: 1 mm/s is w·Δt/Δz = 0.18 at the 3-minute
            # `max_time_step`. On 500 m cells most rivers need one column; the log of `add_rivers`
            # states how many each got.
            maximum_inflow_velocity = 1.0e-3,
            # Water temperature only: every river with a HydAPI daily temperature station on its own
            # watercourse, near the mouth, recorded on all 365 days of 2023. The other mouths carry
            # these rivers' daily mean (`nve_fill_temperatures!`); with none, `add_rivers` fails.
            # Surna's is the one above Trollheim kraftstasjon: the plant's outlet station (112.33.0)
            # lacks 18 days of 2023. Discharge is still every river's REGINE normal; attach
            # `discharge_station`s to the large rivers (Driva, Rauma, Surna, Eira, Aura) as
            # `examples/lista.jl` does.
            outlets = [
                NVERiver(vassdragsnr = "109.A0", name = "Driva", temperature_station = "109.44.0"),
                NVERiver(vassdragsnr = "112.A0", name = "Surna", temperature_station = "112.35.0"),
                NVERiver(vassdragsnr = "103.A20", name = "Rauma", temperature_station = "103.42.0"),
                NVERiver(vassdragsnr = "104.A", name = "Eira", temperature_station = "104.27.0"),
                NVERiver(vassdragsnr = "104.2A0", name = "Visa", temperature_station = "104.23.0"),
                NVERiver(vassdragsnr = "105.A", name = "Osenvassdraget", temperature_station = "105.1.0"),
                NVERiver(vassdragsnr = "096.1A", name = "Hareidselva", temperature_station = "96.3.0"),
                NVERiver(vassdragsnr = "094.4A", name = "Austefjordvassdraget", temperature_station = "94.13.0"),
                NVERiver(vassdragsnr = "102.11Z", name = "Hildreelva", temperature_station = "102.1.0"),
                NVERiver(vassdragsnr = "107.61Z", name = "Sagelva", temperature_station = "107.8.0"),
            ],
        ),
        river_concentrations(),
    ),
    boundary_data_config = EcosmoBoundaries(
        NorKystBoundaryDataConfig(
            data_root        = data_root,
            output_directory = "norkyst_hourly",
            output_file      = "boundaries.nc",
            plot_file        = "boundaries.png",
            # The Norwegian Coastal Current enters through the western edge and leaves through the
            # northern one. In `bathymetry.nc` the western edge is 86 % water and the northern one
            # 85 %; the southern and eastern edges are land but for one fjord-head cell each, 7 m
            # and 5 m deep. The bounding box of the two edge bands is the whole domain, so the
            # hourly download is the whole box for a year: roughly six times Lista's.
            open_edges       = (:west, :north),
            margin           = 0.05,
            architecture     = :auto,
            parameters       = [
                "temperature", "salinity", "u_eastward", "v_northward", "zeta", "ubar", "vbar",
            ],
            years            = years,
        ),
    ),
    atmosphere_config = NORA3Config(
        data_root        = data_root,
        output_directory = "nora3",
        output_file      = "atmosphere.nc",
        plot_file        = "atmosphere.png",
        resolution       = 0.02,
        padding          = 0.1,
        years            = years,
    ),
    model_config = CoupledHydrostaticModel(
        buoyancy           = SeawaterBuoyancy(FT, equation_of_state = TEOS10EquationOfState(FT)),
        closure            = BoundarySponge(
            base = (
                # 7e-7 rather than `oslofjorden.jl`'s 7e-6, as in `oslofjorden_validation.jl`,
                # where 7e-6 mixed the columns towards uniform and freshened the sill basins.
                # Over 640 days at 7e-7 the salinity RMSE against 216 casts was 1.73 g/kg and
                # the temperature RMSE 1.22 °C. Unmeasured here: check `e` against the floor.
                CATKEVerticalDiffusivity(FT; minimum_tke = 7e-7),
                # The validation run's ν = 2e4 m⁴ s⁻¹ scaled by (500/193)⁴ = 45, which keeps its
                # damping times: about 72 min for the 2Δx mode and 56 h at 8Δx. The explicit
                # limit Δt ≤ Δx⁴/32ν₄ is 2000 s at the narrowest 488 m cells.
                #
                # κ is 0, as there: κ = 2e3 left the deep-basin salt budgets unchanged over
                # 60 days but overshot at fronts, leaving 169 wet cells below 0 g/kg.
                HorizontalScalarBiharmonicDiffusivity(FT; ν = 9e5, κ = 0),
            ),
            # The validation run's 16 cells, 8 km here against its 3.2 km, which is under 5 % of
            # this box on each side.
            width_cells = 16,
            # Its 30 and 15 m² s⁻¹ scaled by (500/193)² = 6.7. Δt ≤ Δx²/4ν is 300 s at the
            # narrowest cells, above `max_time_step`.
            viscosity   = 200.0,
            diffusivity = 100.0,
        ),
        # Salinity bounded, as in `oslofjorden_validation.jl`, where plain WENO drove salinity
        # below TEOS10's -32 psu floor beside a large river entering shallow columns. A
        # NamedTuple leaves any tracer it omits on `Centered()`, so it names CATKE's `e` and
        # every ECOSMO tracer as well. ECOSMO's state repair clips their undershoots.
        #
        # Bounded for the three that sink (KNOWN_ISSUES #6): plain WENO drained the cell above the
        # seabed below zero, the state repair clipped it back, and the mass that created piled up in
        # the bottom cell. On this grid's levels at 674.7 m, Float32, 1-minute steps, the bottom
        # cell's detritus went from 5 to 4.2e5 mgC m⁻³ in 7 days with plain WENO, and the run's
        # bottom cells were NaN by day 21; bounded, it settles at 19 with phosphorus conserved to 2e-6.
        # The upper bound is far above any concentration these reach.
        tracer_advection   = merge(
            (T = WENO(FT), S = WENO(FT; bounds = (0, 40)), e = WENO(FT)),
            NamedTuple{ECOSMO_TRACERS}(map(_ -> WENO(FT), ECOSMO_TRACERS)),
            NamedTuple{NivaEcosmo.ECOSMO_SINKING_TRACERS}(
                map(_ -> WENO(FT; bounds = FT.((0, 1e6))), NivaEcosmo.ECOSMO_SINKING_TRACERS),
            ),
        ),
        momentum_advection = WENOVectorInvariant(FT),
        # The whole list: Oceananigans adds a biogeochemistry's tracers only when none of them
        # is named, and `model_tracers` decides the river sources, the open-boundary
        # conditions and the initial conditions before the model exists.
        tracers            = (:T, :S, ECOSMO_TRACERS...),
        coriolis           = HydrostaticSphericalCoriolis(FT),
        sea_ice            = FreezingLimitedOceanTemperature(FT),
        free_surface       = SplitExplicitFreeSurfaceConfig(cfl = 0.7),
        extra_kwargs       = (;),
    ),
    biogeochemistry_config = EcosmoConfig(
        configuration           = ecosmo_configuration,
        surface_shortwave       = surface_shortwave,
        bottom_drag_coefficient = BOTTOM_DRAG_COEFFICIENT,
        start_day_of_year       = START_DAY_OF_YEAR,
    ),
    boundary_condition_configs = (
        AirSeaFluxes(),
        QuadraticBottomDrag(coefficient = BOTTOM_DRAG_COEFFICIENT),
        # Marchesiello et al. (2001): nudge hard on inflow, radiate on outflow. Oslofjord's
        # timescales, applied to every open edge.
        OpenLateralBoundaryFromData(
            inflow_timescale  = 3hours,
            outflow_timescale = 360days,
        ),
        EcosmoGasExchange(
            wind_speed       = WIND_SPEED,
            atmospheric_pCO2 = ATMOSPHERIC_PCO2,
            pH_scale         = ecosmo_configuration.pH_scale,
        ),
    ),
    run_config = RunConfig(
        results_root       = fjord_results_root("moro"),
        architecture       = :auto,
        writers = (
            # Daily, as in the validation run. One record of these five fields is 147 MB, so the
            # year is 54 GB.
            SnapshotWriter(
                name               = :ocean,
                output_file        = "snapshots_ocean.nc",
                variables          = (:T, :S, :u, :v, :e),
                interval           = 1day,
                overwrite_existing = true,
            ),
            # Weekly, for size: one record of the 17 tracers is 500 MB, so 26 GB over the year
            # against 180 GB daily.
            SnapshotWriter(
                name               = :ecosmo,
                output_file        = "snapshots_ecosmo.nc",
                variables          = ECOSMO_TRACERS,
                interval           = 7days,
                overwrite_existing = true,
            ),
            # Hourly, one file per calendar month. 35 surface fields of 396 x 315 Float32 values are
            # 17.5 MB a record before compression, 153 GB over the year; `JLD2Writer` deflates by
            # default.
            EcosmoSurfaceWriter(
                name               = :surface,
                output_file        = "snapshots_surface.jld2",
                interval           = 1hour,
                split_times        = MONTH_STARTS,
                overwrite_existing = true,
            ),
            # The two permanent Kartverket gauges in the box, from Kartverket's station register
            # (`vannstand.kartverket.no/tideapi.php?tide_request=stationlist&type=perm`). Ten
            # minutes, as in the validation run, for the harmonic analysis `validate_simulation`
            # does against Kartverket's own record.
            FieldStationWriter(
                name               = :tidegauge,
                output_file        = "stations_tidegauge.jld2",
                variables          = (:η,),
                stations           = [
                    Station(name = "Ålesund",      longitude = 6.151946, latitude = 62.469414),
                    Station(name = "Kristiansund", longitude = 7.736140, latitude = 63.113920),
                ],
                interval           = 10minutes,
                overwrite_existing = true,
            ),
            CheckpointWriter(interval = 12hours, cleanup = true),
        ),
        callbacks = (ProgressCallback(name = :progress, interval = 1hour, report = progress),),
        # `oslofjorden.jl`'s. The advective CFL limits Δt, not `max_time_step`: the validation run
        # took 25-32 s steps at 193 m, so this grid should take roughly 2.6 times that.
        time_stepping = AdaptiveTimeStep(
            initial_time_step    = 1second,
            cfl                  = 0.3,
            max_time_step        = 3minutes,
            max_time_step_change = 1.01,
        ),
        # The NorKyst daily mean of `start_date` for T, S, u and v, and the ECOSMO profiles.
        initial_conditions = WithProfiles(FromInitialState(), ECOSMO_PROFILES),
        # The calendar year 2023.
        start_date         = START_DATE,
        stop_time          = 365days,
        loops              = 1,
        pickup             = false,
    ),
)
