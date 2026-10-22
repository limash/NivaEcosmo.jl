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
| `src/parameters.jl` | `EcosmoParameters`: every parameter with its Fortran name, default and FABM unit conversion |
| `src/ecosmo.jl` | `interior_rates` (`do`), `oxygen_surface_flux` (`do_surface`), `bottom_fluxes` (`do_bottom`), `sinking_speeds` |
| `src/carbonate.jl` | `co2_constants` (`CO2SET`), `carbonate_system` (`CO2DYN`/`CO2CLC`), `calcite_aragonite_saturation`, the carbonate `do` and `do_surface` |
| `src/daylength.jl`, `src/light.jl` | day length; `gotm/light` and the `GOTMLight` OceanBioME component |
| `src/biogeochemistry.jl` | `ECOSMO`, the OceanBioME model: tracer tendencies, sinking, auxiliary fields, `ClipToMinimum` |
| `src/sediment.jl` | `ecosmo_sediment`: `sed1`–`sed4` as an OceanBioME `BiogeochemicalSediment` |
| `src/boundary_conditions.jl` | `ecosmo_surface_fluxes`: O₂ and CO₂ air-sea exchange as flux boundary conditions |
| `src/fabm_yaml.jl` | `fabm_configuration`: parameters, light and initial values from a niva-ecosmo `fabm.yaml` |
| `box_model/` | the Drammensfjorden box model and its comparison with the Fortran (see its README) |
| `test/runtests.jl` | tests |

The processes are pointwise functions of one cell, shared by the OceanBioME kernels and the box
model. They keep the Fortran's evaluation order, so the two implementations round alike.

## Does it match the Fortran?

- **Box model, one year** (`box_model/`, the twin of niva-ecosmo's `run_box.py`): every state
  variable and diagnostic agrees with pyfabm to within 3.3 × 10⁻¹⁵ of its range. Given the same
  state, rates agree to a few units in the last place, and most are bit-identical. The remainder is
  `exp`/`log`/`pow`/trigonometry rounding: Julia implements those itself, gfortran takes them from
  glibc.
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

configuration = fabm_configuration("box_model/fabm.yaml")   # or EcosmoParameters(; muPs = 1.47, ...)
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
                          ecosmo_surface_fluxes(; wind_speed = 5.0, pH_scale = configuration.pH_scale))

model = HydrostaticFreeSurfaceModel(grid; tracers = (:T, :S, ECOSMO_TRACERS...), biogeochemistry,
                                    boundary_conditions, buoyancy = SeawaterBuoyancy())

set!(model; configuration.initial_values[ECOSMO_TRACERS]..., T = 8, S = 30)
for name in (:sed1, :sed2, :sed3, :sed4)
    set!(biogeochemistry.sediment.fields[name], configuration.initial_values[name])
end
```

The tracers keep their FABM names and units (`ECOSMO_TRACERS`): `no3 nh4 pho sil oxy fla dia
flachl diachl microzoo mesozoo det opa cal dom dic alk`, nutrients, plankton and organic matter in
mgC m⁻³ (Redfield), chlorophyll in mg m⁻³, `oxy`, `dic` and `alk` in mmol m⁻³. The model reads `T`
(°C) and `S` (PSU) but does not require them as biogeochemical tracers. The auxiliary fields `PAR`,
`pH`, `pCO2`, `Carb`, `Om_cal` and `day_length` can be written like any model field.

### From FjordSim

The same three steps as `examples/oslofjorden_npzd.jl`, with `ECOSMO` in place of `NPZD`:

```julia
using NivaEcosmo
using Dates: dayofyear

struct EcosmoModel{M,T,S,W} <: AbstractCoupledSimulationConfig
    base::M
    tracers::T                  # (:T, :S, ECOSMO_TRACERS...): the list must be complete
    surface_shortwave::S
    wind_speed::W
end

Simulations.model_tracers(model::EcosmoModel) = model.tracers

Simulations.coupled_simulation(model::EcosmoModel, grid; kwargs...) = Simulations.coupled_simulation(
    CoupledHydrostaticSimulation(
        # ... every other field from model.base, as npzd_model does ...
        tracers = model.tracers,
        biogeochemistry = ECOSMO(grid;
                                 surface_shortwave = model.surface_shortwave,
                                 start_day_of_year = dayofyear(start_date) - 1,
                                 sediment = ecosmo_sediment(grid; bottom_stress = QuadraticBottomStress(drag_coefficient = 0.003))),
    ),
    grid; kwargs...)

# the gas exchange, as one more piece of MergedBoundaryConditions
struct EcosmoSurfaceFluxes{W} <: AbstractBoundaryConditionConfig
    wind_speed::W
end

BoundaryConditions.boundary_condition_sides(config::EcosmoSurfaceFluxes, grid, forcing, boundary_config, tracers) =
    map(bc -> (top = bc,), ecosmo_surface_fluxes(; wind_speed = config.wind_speed))
```

As with NPZD, FjordSim's open boundaries need a prepared series for every one of the 17 tracers,
and the rivers may carry any of them.

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
