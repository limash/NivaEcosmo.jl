# NivaEcosmo.jl

The NIVA adaptation of the ECOSMO marine biogeochemical model
([niva-ecosmo](../niva-ecosmo), FABM/Fortran) ported to Julia as an
[OceanBioME](https://github.com/OceanBioME/OceanBioME.jl) biogeochemistry, for
[Oceananigans](https://github.com/CliMA/Oceananigans.jl) models and so for
[FjordSim](../FjordSim.jl).

Ported from niva-ecosmo commit `b351e6f`: `niva_ecosmo.F90` (pelagic and benthic processes),
`niva_ecosmo_carbonate.F90` (carbonate system and air-sea CO₂), `niva_ecosmo_daylength.F90`, and
FABM's `gotm/light`, the light model the configurations run with.

## Layout

| Path | Content |
|---|---|
| `src/constants.jl` | `ECOSMO_CONSTANTS` and `ecosmo_constants(FT)`: every fixed number of the processes, with its source and units |
| `src/parameters.jl` | `EcosmoParameters`: every parameter with its Fortran name, default and FABM unit conversion |
| `src/ecosmo.jl` | `interior_rates` (`do`), `oxygen_surface_flux` (`do_surface`), `bottom_fluxes` (`do_bottom`), `sinking_speeds` |
| `src/carbonate.jl` | `co2_constants` (`CO2SET`), `carbonate_system` (`CO2DYN`/`CO2CLC`), `calcite_aragonite_saturation`, the carbonate `do` and `do_surface` |
| `src/daylength.jl`, `src/light.jl` | day length; `gotm/light` and the `GOTMLight` OceanBioME component |
| `src/biogeochemistry.jl` | `ECOSMO`, the OceanBioME model: tracer tendencies, sinking, auxiliary fields, `ClipToMinimum` |
| `src/sediment.jl` | `ecosmo_sediment`: `sed1`–`sed4` as an OceanBioME `BiogeochemicalSediment` |
| `src/boundary_conditions.jl` | `ecosmo_surface_fluxes`: O₂ and CO₂ air-sea exchange as flux boundary conditions |
| `src/fabm_yaml.jl` | `fabm_configuration`: parameters, light and initial values from a niva-ecosmo `fabm.yaml` |
| `data/` | input files: `fabm.yaml` (the A20 v3aa configuration the tests and the box model read), the box model's forcing and boundary values, and the FjordSim setups' tracer values |
| `box_model/` | the Drammensfjorden box model and its comparison with the Fortran (see its README) |
| `fjordsim_model/` | FjordSim setups running ECOSMO, with their own environment (see its README) |
| `test/runtests.jl` | tests |

The processes are pointwise functions of one cell, shared by the OceanBioME kernels and the box
model. Every number in them is a parameter, an exact integer or an entry of `ECOSMO_CONSTANTS`
(`src/constants.jl`), so they compute in the float type of the parameters and constants they are
given.

## Does it match the Fortran?

- **Box model, one year** (`box_model/`, the twin of niva-ecosmo's `run_box.py`): every state
  variable and diagnostic agrees with pyfabm to within 1.1 × 10⁻¹⁴ of its range, the carbonate
  diagnostics being the furthest. The two round differently in the last places: the formulas are
  not written in the Fortran's form, and Julia implements `exp`, `log`, `pow` and the trigonometric
  functions itself where gfortran takes them from glibc.
- **The OceanBioME model**: in a one-cell Oceananigans model, every tracer tendency equals the box
  model's interior rate plus surface and bottom flux over the cell thickness (relative 10⁻¹²), and
  the sediment tendencies are bit-identical (`test/runtests.jl`).
- **3D**: on a latitude-longitude grid with an immersed seabed, with sinking, gas exchange,
  sediment and resuspension, water-column plus sediment phosphorus is conserved to round-off
  (`test/runtests.jl`).

```bash
julia --project test/runtests.jl
```

## Use in an Oceananigans model

```julia
using Oceananigans, NivaEcosmo

configuration = fabm_configuration("data/fabm.yaml", eltype(grid))  # or EcosmoParameters(eltype(grid); ...)
parameters = configuration.parameters

biogeochemistry = ECOSMO(grid;
                         parameters,
                         surface_shortwave = (λ, φ, t) -> 200.0,  # W/m², shortwave not PAR
                         light_attenuation = GOTMLight(grid, (λ, φ, t) -> 200.0; configuration.light...),
                         sediment = ecosmo_sediment(grid; parameters,
                                                    bottom_stress = QuadraticBottomStress(drag_coefficient = 0.003)),
                         start_day_of_year = 0,                   # model time zero is 1 January 00:00
                         pH_scale = configuration.pH_scale)

boundary_conditions = map(bc -> FieldBoundaryConditions(top = bc),
                          ecosmo_surface_fluxes(eltype(grid); wind_speed = 5.0, pH_scale = configuration.pH_scale))

model = HydrostaticFreeSurfaceModel(grid; tracers = (:T, :S, ECOSMO_TRACERS...), biogeochemistry,
                                    boundary_conditions, buoyancy = SeawaterBuoyancy())

set!(model; configuration.initial_values[ECOSMO_TRACERS]..., T = 8, S = 30)
for name in (:sed1, :sed2, :sed3, :sed4)
    set!(biogeochemistry.sediment.fields[name], configuration.initial_values[name])
end
```

The float type is the grid's: `ECOSMO` takes it from the grid, so a `Float32` grid needs
`fabm_configuration(path, Float32)` or `EcosmoParameters(Float32; ...)`, and
`ecosmo_surface_fluxes(Float32; ...)`.

### Float32

On a `Float32` grid the whole biogeochemistry computes in `Float32`: the parameters and
`ecosmo_constants(Float32)` are converted once, when the model is built. A shortwave or wind
function that returns `Float64` is converted where a kernel reads it.

Against `Float64`, measured on 2000 random states and in `test/runtests.jl`'s 3D model:

- Rates and bottom fluxes agree to within 7 × 10⁻⁷ of the cell's largest rate. After twenty
  10-minute steps of the 3D model every tracer and sediment pool is within 7 × 10⁻⁷ of its largest
  value, and total phosphorus is conserved to 3 × 10⁻⁸.
- The carbonate system loses more. Its equilibrium constants are logarithms of polynomials whose
  terms are up to a thousand times their sum, so the cancellation magnifies single-precision
  rounding: `pCO2`, `Carb` and `Om_cal` are off by up to 2 × 10⁻⁴ relative (0.08 µatm at
  400 µatm), `pH` by 10⁻⁵. The iteration converged in all 6000 states tried, 2000 on each pH scale.
- The air-sea fluxes are differences between air and water: the CO₂ flux error is that of the
  pCO₂ error, the O₂ error that of a 0.007 mmol m⁻³ departure from saturation.

The tracers keep their FABM names and units (`ECOSMO_TRACERS`): `no3 nh4 pho sil oxy fla dia
flachl diachl microzoo mesozoo det opa cal dom dic alk`, nutrients, plankton and organic matter in
mgC m⁻³ (Redfield), chlorophyll in mg m⁻³, `oxy`, `dic` and `alk` in mmol m⁻³. The model reads `T`
(°C) and `S` (PSU) but does not require them as biogeochemical tracers. The auxiliary fields `PAR`,
`pH`, `pCO2`, `Carb`, `Om_cal` and `day_length` can be written like any model field.

### From FjordSim

`fjordsim_model/moro.jl` is a complete FjordSim setup running ECOSMO; `fjordsim_model/README.md`
says how to install and run it. It follows FjordSim's `examples/oslofjorden_npzd.jl`, through the
hooks FjordSim offers an out-of-tree setup:

- **The model.** An `AbstractBiogeochemistryConfig`, the setup's `biogeochemistry_config`, whose
  `Simulations.model_biogeochemistry(config, grid)` method builds `ECOSMO`, because `ECOSMO` needs
  the grid and a setup file is evaluated before there is one. The sediment's pools are set there too, since they are not model fields.
  `tracers` must name all of them, `(:T, :S, ECOSMO_TRACERS...)`, and so must a `NamedTuple`
  `tracer_advection`: a tracer it omits is advected by `Centered()`.
- **Gas exchange.** An `AbstractBoundaryConditionConfig` whose `boundary_condition_sides` returns
  `ecosmo_surface_fluxes` as top conditions, one more piece of `boundary_condition_configs`.
  FjordSim's `AirSeaFluxes` gives a top condition to T and S only.
- **Open boundaries.** `OpenLateralBoundaryFromData` radiates every model tracer towards its
  prepared series and raises for a tracer with none, so an `AbstractBoundaryDataConfig` wrapping the
  NorKyst one appends all 17 to the prepared boundary file.
- **Rivers.** A tracer the river file does not carry enters at zero. An `AbstractRiverConfig`
  wrapping `NVERiversConfig` adds a series per tracer to its `river_series`.
- **Initial state.** `FromInitialState` reads T, S, u and v only, so a wrapper of it adds the
  ECOSMO tracers through `Simulations.resolve_initial_conditions`.

## Where the port differs from the Fortran

What the Fortran gets from its FABM host has to come from somewhere in Oceananigans. In the box
model all of it is given exactly as pyfabm gives it; in a 3D model:

- **Density** for the carbonate unit conversions is `density(T, S, P)`, by default OceanBioME's
  TEOS-10 polynomial; **pressure** is the depth in metres, taken as dbar.
- **Bottom stress** is a constant (default 0: everything settles, nothing is resuspended) or
  `QuadraticBottomStress`, ρ₀ C_d |u|² from the bottom-cell velocity.
- **Day of year** runs on a 365-day calendar from `start_day_of_year`, and the latitude comes from
  the grid or `latitude`.
- **Wind and atmospheric pCO₂** for the gas exchange are inputs of `ecosmo_surface_fluxes`.
- **State repair**: FABM hosts clip state variables to their minimum (`repair_state`). `ECOSMO`
  does the same before each tendency evaluation with `ClipToMinimum`, on by default, because the
  chlorophyll equations divide by phytoplankton biomass. It is not mass conserving.
- **Sinking** is advective in the interior, at `sinking_speeds` for each face's depth, and zero on
  the seabed and the surface, an immersed seabed included. Material leaves the bottom cell by the
  sediment's settling flux, as in FABM. That flux uses the minimum speed `sinkD`, not the
  depth-dependent one, also as in the Fortran.
- **Sediment fluxes** are divided by the bottom cell's thickness. OceanBioME's default divides by
  the distance between cell centres, which does not conserve mass on a stretched grid or a partial
  cell, so `ecosmo_sediment` overrides it.

Not ported: `use_aice` (ice-limited air-sea fluxes) and `couple_co2: false`. `fabm_configuration`
rejects both. Also not ported are the `TSOM_s`/`TSIM_s` suspended-matter aggregates, which only
`light_iop_aksnes.F90` reads, and that is not part of the niva-ecosmo build.

Issues found in the Fortran and in the Julia packages, and a decision on light attenuation that is
open, are in [KNOWN_ISSUES.md](KNOWN_ISSUES.md).

## Performance

Each tracer's tendency evaluates the whole interior model and keeps its own rate, which is how
OceanBioME's models are written. What a tracer does not use is only partly compiled away, so this
is the first thing to revisit if the biogeochemistry dominates a run's cost. The pointwise
functions are type-stable and allocation-free. No GPU was available to test on.

On GPUs whose `Float64` throughput is a small fraction of `Float32` (1/32 to 1/64 on most
consumer and workstation cards, 1/2 on data-centre ones), run on a `Float32` grid: ECOSMO's
arithmetic then stays in `Float32`. On a CPU, Julia's own `Float32` `^`, `log`, `log10` and `cos`
compute internally in `Float64`; CUDA.jl replaces them with single-precision device functions.
