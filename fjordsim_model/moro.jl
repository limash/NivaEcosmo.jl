# Møre og Romsdal with ECOSMO: FjordSim's `examples/moro/moro.jl` (the coast from Stad to Hitra,
# 5.12-9.00°E by 62.00-63.43°N at ~500 m, 52 levels to 1137 m, NorKyst-800m boundaries on the western
# and northern edges, NVE rivers, NORA3, the year 2023), with ECOSMO's 17 tracers beside T and S and
# the A20 v3aa parameters of `data/fabm.yaml`. The FjordSim parts are NivaEcosmo's FjordSim extension
# (`ext/NivaEcosmoFjordSimExt.jl`); this file states the values.
#
#   julia --project=fjordsim_model -m FjordSim <step> --config fjordsim_model/moro.jl
#
# Installation and the order of the steps are in `fjordsim_model/README.md`. Data goes to
# `~/FjordSim_data/moro/`, shared with the FjordSim setup, whose downloads, bathymetry and atmosphere
# serve this file too; the initial state, the rivers and the boundaries carry the tracers, so this
# file prepares its own beside them (`*_ecosmo.nc`). Results go to `~/FjordSim_results/moro/`.
#
# **Every biogeochemical value is a placeholder, not a measurement.** They are in
# `data/moro_ecosmo.yaml`: typical values for Atlantic-influenced Norwegian coastal water and
# west-coast rivers, with a seasonal cycle, uniform along the edges. Replace them with observations
# or a reanalysis (FjordSim's `examples/roho/roho_ecosmo.jl` reads TOPAZ5-ECOSMO) before the run is
# used for anything but bringing the setup up. Light is a clear-sky shortwave scaled by a constant
# clearness index, and the gas exchange runs at a constant wind.

using FjordSim
using NivaEcosmo
using NivaEcosmo: ECOSMO_SINKING_TRACERS
# Loads NivaEcosmo's FjordSim extension, which writes the tracers into the prepared files.
using NCDatasets: NCDatasets
using Oceananigans: WENO
using Oceananigans.Units: day, days, hour
using Dates: DateTime, Second, dayofyear
using YAML

base = include(joinpath(pkgdir(FjordSim), "examples", "moro", "moro.jl"))

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
logistic opening at `onset_day` times one closing at `end_day`.
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

Every ECOSMO tracer as the extension's `(z, date, source) -> value`, ignoring `source`. A tracer the
file does not give is a `KeyError` here, rather than a missing boundary series at the start of a run.
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
    return NamedTuple{ECOSMO_TRACERS}(map(name -> (z, date, source) -> profiles[name](z, date), ECOSMO_TRACERS))
end

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

# --- Light, wind and CO₂ ------------------------------------------------------------------------

# Days since midnight on 1 January, as ECOSMO's `start_day_of_year` counts them.
const START_DAY_OF_YEAR = dayofyear(start_date) - 1
# The clear-sky shortwave of `examples/oslofjorden/oslofjorden_npzd.jl`, as shortwave rather than
# PAR: `GOTMLight` keeps its own visible fraction. The clearness index is that file's value for Oslo,
# so the cloudier west coast likely gets too much light.
const SOLAR_CONSTANT = 1361.0
const CLEARNESS_INDEX = 0.5
# 10 m wind for the O₂ and CO₂ exchange, m s⁻¹, assumed. The exchange goes as the square of the
# wind, so a mean speed underestimates it.
const WIND_SPEED = 7.0
# CO₂ in air, ppm: the 2023 global mean to within about 1 ppm.
const ATMOSPHERIC_PCO2 = 420.0

"""
    surface_shortwave(λ, φ, t)

Daily-mean downwelling shortwave below the surface, W m⁻², at `t` seconds after `start_date`: the
top-of-atmosphere daily mean (Cooper's declination) times `CLEARNESS_INDEX`. Runs in a kernel.
"""
@inline function surface_shortwave(λ, φ, t)
    day_of_year = START_DAY_OF_YEAR + t / 86400
    declination = deg2rad(23.44) * sin(2π * (day_of_year + 284) / 365)
    latitude = deg2rad(φ)
    hour_angle = acos(clamp(-tan(latitude) * tan(declination), -1, 1))
    eccentricity = 1 + 0.033 * cos(2π * day_of_year / 365)
    insolation = SOLAR_CONSTANT / π * eccentricity *
                 (hour_angle * sin(latitude) * sin(declination) + cos(latitude) * cos(declination) * sin(hour_angle))
    return CLEARNESS_INDEX * max(0, insolation)
end

# --- The setup ----------------------------------------------------------------------------------

# The A20 v3aa parameters, unchanged, with the carbonate module's pH scale and the sediment's
# initial pools, in the grid's `FT`, which `ECOSMO` requires.
ecosmo_configuration = fabm_configuration(joinpath(pkgdir(NivaEcosmo), "data", "fabm.yaml"), FT)

const PROFILES = ecosmo_profiles()

# The FjordSim setup's components, writing prepared files of their own.
function renamed(config; files...)
    renamed_config = deepcopy(config)
    for (field, name) in pairs(files)
        setproperty!(renamed_config, field, name)
    end
    return renamed_config
end

physics = base.model_config

FjordConfig(
    grid_config = base.grid_config,
    bathymetry_config = base.bathymetry_config,
    initial_state_config = ecosmo_initial_state(
        renamed(base.initial_state_config; output_file = "initial_state_ecosmo.nc");
        tracers = PROFILES,
    ),
    river_config = ecosmo_rivers(
        renamed(base.river_config; output_file = "rivers_ecosmo.nc", plot_file = "rivers_ecosmo.png"),
        river_concentrations(),
    ),
    boundary_data_config = ecosmo_boundaries(
        renamed(base.boundary_data_config; output_file = "boundaries_ecosmo.nc", plot_file = "boundaries_ecosmo.png");
        tracers = PROFILES,
    ),
    atmosphere_config = base.atmosphere_config,
    model_config = CoupledHydrostaticModel(
        buoyancy = physics.buoyancy,
        closure = physics.closure,
        # The ECOSMO tracers on plain WENO, the state repair clipping their undershoots, but bounded
        # for the three that sink (KNOWN_ISSUES #6): plain WENO drained the cell above the seabed
        # below zero, the repair clipped it back, and the mass that created piled up in the bottom
        # cell. On this grid's levels at 674.7 m, Float32, 1-minute steps, the bottom cell's detritus
        # went from 5 to 4.2e5 mgC m⁻³ in 7 days with plain WENO, and the run's bottom cells were NaN
        # by day 21; bounded, it settles at 19 with phosphorus conserved to 2e-6.
        tracer_advection = merge(
            physics.tracer_advection,
            NamedTuple{ECOSMO_TRACERS}(map(_ -> WENO(FT), ECOSMO_TRACERS)),
            NamedTuple{ECOSMO_SINKING_TRACERS}(map(_ -> WENO(FT; bounds = FT.((0, 1e6))), ECOSMO_SINKING_TRACERS)),
        ),
        momentum_advection = physics.momentum_advection,
        # The whole list: the river sources, the open-boundary conditions and the initial conditions
        # are decided from it before the model exists.
        tracers = (:T, :S, ECOSMO_TRACERS...),
        coriolis = physics.coriolis,
        sea_ice = physics.sea_ice,
        free_surface = physics.free_surface,
        extra_kwargs = physics.extra_kwargs,
    ),
    biogeochemistry_config = ecosmo_biogeochemistry_config(
        configuration = ecosmo_configuration,
        surface_shortwave = surface_shortwave,
        # The FjordSim setup's bottom drag, which the sediment's resuspension stress should share.
        bottom_drag_coefficient = 0.003,
        start_day_of_year = START_DAY_OF_YEAR,
    ),
    boundary_condition_configs = (
        base.boundary_condition_configs...,
        ecosmo_gas_exchange(
            wind_speed = WIND_SPEED,
            atmospheric_pCO2 = ATMOSPHERIC_PCO2,
            pH_scale = ecosmo_configuration.pH_scale,
        ),
    ),
    run_config = RunConfig(
        results_root = fjord_results_root("moro"),
        architecture = base.run_config.architecture,
        writers = (
            # The FjordSim setup's writers but its surface one, which the ECOSMO surface writer
            # replaces.
            filter(writer -> !(writer isa SurfaceSnapshotWriter), base.run_config.writers)...,
            # Weekly, for size: one record of the 17 tracers is 500 MB, so 26 GB over the year
            # against 180 GB daily.
            SnapshotWriter(
                name = :ecosmo,
                output_file = "snapshots_ecosmo.nc",
                variables = ECOSMO_TRACERS,
                interval = 7days,
                split_times = Float64[],
                overwrite_existing = true,
            ),
            # Hourly, one file per calendar month. 35 surface fields of 396 x 315 Float32 values are
            # 17.5 MB a record before compression, 153 GB over the year; `JLD2Writer` deflates.
            ecosmo_surface_writer(
                name = :surface,
                output_file = "snapshots_surface.jld2",
                interval = 1hour,
                split_times = [(DateTime(2023, month, 1) - DateTime(2023, 1, 1)) / Second(1) for month in 2:12],
                overwrite_existing = true,
            ),
        ),
        callbacks = base.run_config.callbacks,
        time_stepping = base.run_config.time_stepping,
        # The NorKyst daily mean of `start_date` for T, S, u and v, and the ECOSMO profiles.
        initial_conditions = FromInitialState(),
        start_date = start_date,
        stop_date = stop_date,
        loops = 1,
        pickup = false,
    ),
)
