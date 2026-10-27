# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Docs — read before the matching task

| Doc | Read it before |
|---|---|
| [`README.md`](README.md) | using `ECOSMO` in a model or from FjordSim; anything touching where the port differs from the Fortran |
| [`KNOWN_ISSUES.md`](KNOWN_ISSUES.md) | changing light attenuation, settling/sinking, or the sediment coupling |
| [`box_model/README.md`](box_model/README.md) | running the box model or comparing against pyfabm |
| [`fjordsim_model/README.md`](fjordsim_model/README.md) | installing or running the FjordSim setups (`moro.jl`) |

## Commands

```bash
julia --project test/runtests.jl                                        # all tests
julia --project=box_model -e 'using Pkg; Pkg.instantiate()'             # box model env, once
julia --project=box_model box_model/run_box.jl                          # one year -> box_model/output.nc, output.png
julia --project=box_model box_model/run_box.jl --closed --no-plot       # no exchange, no burial: total P conserved
julia --project=box_model box_model/compare.jl <pyfabm output.nc>       # compare with niva-ecosmo's box_model output
julia --project=fjordsim_model -e 'using Pkg; Pkg.instantiate()'        # FjordSim setups env, once
julia --project=fjordsim_model -m FjordSim <step> --config fjordsim_model/moro.jl
```

Input files — `fabm.yaml`, the box model's `forcing.csv` and `boundary.yaml`, and the FjordSim
setups' tracer values (`moro_ecosmo.yaml`) — are in `data/`.

The tests are one file of nested `@testset`s with no filter; to run one, evaluate the file's
`using` block and that `@testset` in a REPL started with `julia --project`. `.JuliaFormatter.toml`
sets a 120-character margin.

## Architecture

A port of niva-ecosmo (FABM/Fortran, commit `b351e6f`) to an OceanBioME biogeochemistry. Two
layers:

1. **Pointwise process functions** of one cell, one per Fortran routine: `interior_rates` (`do`),
   `bottom_fluxes` (`do_bottom`), `oxygen_surface_flux` (`do_surface`), `sinking_speeds`
   (`src/ecosmo.jl`); `carbonate_diagnostics`, `carbon_dioxide_surface_flux` (`src/carbonate.jl`);
   `day_length`; `gotm_light`. They take the parameters and `ecosmo_constants(FT)`
   (`src/constants.jl`) and return NamedTuples (`(; rates, diagnostics)`, `(; fluxes, sediment)`)
   keyed by FABM names.
2. **Two drivers** call them: the OceanBioME wiring in `src/biogeochemistry.jl`, `src/sediment.jl`,
   `src/light.jl` and `src/boundary_conditions.jl`, and the box model (`box_model/run_box.jl`),
   which repeats pyfabm's `run_box.py` step for step. A process change goes in the pointwise
   function only, so both drivers see it.

How the OceanBioME side is put together:

- `ECOSMO(grid; ...)` builds `EcosmoBiogeochemistry` and wraps it in OceanBioME's `Biogeochemistry`
  with `GOTMLight`, `ecosmo_sediment` and the `ClipToMinimum` modifier (FABM's `repair_state`, on
  by default — the chlorophyll equations divide by phytoplankton carbon).
- Each tracer's tendency method calls `cell_rates`, which evaluates the whole `interior_rates` and
  keeps one field of it (methods generated with `@eval` over `ECOSMO_TRACERS`).
- Diagnostics the rates depend on — day length and the carbonate system (`pH`, `pCO2`, `Carb`,
  `Om_cal`) — are auxiliary fields computed in `update_biogeochemical_state!` before any tendency,
  as FABM does. The carbonate fields carry state between steps: a failed carbonate iteration falls
  back to the previous values.
- Sinking (`det`, `opa`, `cal`) is a drift velocity field, zero on the seabed and surface (immersed
  seabed included); material leaves the bottom cell only through the sediment's settling flux.
- `src/sediment.jl` overrides OceanBioME's `update_tendencies!` for the ECOSMO sediment to divide by
  `Δzᶜᶜᶜ` (KNOWN_ISSUES #5).
- Air-sea O₂ and CO₂ are flux boundary conditions from `ecosmo_surface_fluxes`, attached by the
  user, not by `ECOSMO`.
- `EcosmoParameters` (`src/parameters.jl`) is generated from the `ECOSMO_PARAMETERS` table of
  `(name, default, scale factor, unit, description)`; values are stored already scaled (per
  second, mgC m⁻³), exactly as FABM's `get_parameter(..., scale_factor)` stores them.
  `fabm_configuration` reads a niva-ecosmo `fabm.yaml` into parameters, light coefficients and
  initial values.

## Precision and agreement with the Fortran

- The float type is chosen once, by the grid: `ECOSMO` builds `ecosmo_constants(FT)` from it, and
  `EcosmoParameters(FT)`, `fabm_configuration(path, FT)` and `ecosmo_surface_fluxes(FT; ...)` take
  the same `FT`. The box model is `Float64`.
- A pointwise function contains no decimal literals. Every number in it is a parameter, an entry
  of `ECOSMO_CONSTANTS` (`src/constants.jl`, with its source and units), or an exact integer
  (`1 - x`, `x / 2`, `(temp - 10) / 10`). A new number goes in the table. Integer division `1 / 2`
  is a `Float64`.
- In `Float64` the box model agrees with pyfabm to ~10⁻¹⁴ of each variable's range; the `Float32`
  agreement is in the README ("Float32"). The formulas need not copy the Fortran's form.
- Parameter names are the `fabm.yaml` keys, and `interior_rates` and `bottom_fluxes` keep the
  Fortran's variable names (`muPs`, `up_nh4Ps`, `E0c`) so a process can be found in
  `niva_ecosmo.F90`.
- Behaviour copied from the Fortran on purpose, including KNOWN_ISSUES #1–3 and the oddities noted
  in `ECOSMO_CONSTANTS` (two single-precision coefficients, two compressibilities not divided by
  1000), is not "fixed" in Julia alone: a change there is made in both implementations or the two
  stop matching.
- After touching a process, rerun the tests and the box model comparison.

## GPU & kernel compatibility

Applies to `@kernel` functions and everything they call (`interior_rates`, `bottom_fluxes`, the
carbonate functions, tendency and sediment methods):

- Mark with `@inline`; keep type-stable and allocation-free.
- In new kernel code use `ifelse`, not `&&`/`||`/`?:`. The existing branches in the ported
  functions (`if E0 > constants.minimum_light`, the carbonate iteration) mirror the Fortran's
  control flow — `ifelse` evaluates both sides — and stay as they are.
- Never loop over grid points outside a kernel — use `launch!`.
- Use `on_architecture` for CPU↔GPU transfers — never `Array()`/`CuArray()` directly.
- Index fields with three indices, `field[i, j, k]`; 2D fields as `field[i, j, 1]`.
- A struct that a kernel receives needs `Adapt.adapt_structure` for its field members (see
  `EcosmoBiogeochemistry`).
- Kernel code computes in the grid's `FT`: with the parameters and constants it is given,
  `zero(...)`, `one(...)` and integers, never decimal literals. A user's function that may return
  `Float64` (wind, shortwave) is converted to `eltype(grid)` where the kernel reads it. The tests
  check `Float32` result types and a `Float32` run against a `Float64` one. Nothing has been run
  on a GPU (README, "Performance").

## Type stability

- Structs have concrete field types — parametrize rather than use `Any` or abstract fields.
- User-facing constructors may take flexible arguments but return fully typed structs (as
  `ECOSMO` converts `latitude` and `start_day_of_year` to `FT`, `ecosmo_sediment` the stress).
- Type annotations are for dispatch, not documentation.
- Check `@code_warntype` when touching a kernel-executed path.

## Import conventions

- Exports go at the top of `src/NivaEcosmo.jl`, before anything else. Export only what a user
  names in a model setup; internals (`interior_rates`, `carbonate_system`, constants) are reached
  as `NivaEcosmo.name`, as the tests do.
- Imports are explicit (`using Mod: a, b`), all in `src/NivaEcosmo.jl`, so the dependency surface
  stays auditable. Functions extended from Oceananigans and OceanBioME are listed in its
  `import ...: ...` blocks; add a new extended function there.

## Tests

- Tests are offline and build their own fixtures; parameters and initial values come from
  `data/fabm.yaml`, a file copied unchanged from niva-ecosmo.
- Test a property the port must have — conservation, equality of the OceanBioME tendencies with
  the pointwise composition, Fortran agreement — not an incidental number.
- Do not test plotting.

## Code style

- Match the surrounding code: spaces around `=` in keyword arguments (`EcosmoParameters(; muPs =
  2.0)`), 120-character lines.
- In new Julia-side code, avoid abbreviations (`latitude`, not `lat`; physics symbols `T`, `S`,
  `u`, `v`, `λ`, `φ` are fine). Ported code keeps the Fortran names (see above).
- No trailing whitespace or blank lines; files end with exactly one newline.
- Explicit `return` in functions longer than one expression.
- Delete commented-out code and debugging artifacts.
- Comments state what holds for any reader on any machine: what the code does and why, with the
  Fortran routine or line it ports, and measured numbers with what they were measured on. Never
  the author's situation or the file's history.

## Writing docs and CLAUDE.md

- No user- or account-specific information in any checked-in file: no names, email addresses, one
  machine's local paths, or one person's preferences. Those go in Claude's memory.
- CLAUDE.md holds only what every task needs; rationale goes in `README.md` or `KNOWN_ISSUES.md`.
- Keep `KNOWN_ISSUES.md` current: when an issue is resolved or decided, update or remove its entry.

## Common pitfalls

- Never set `Oceananigans.defaults.FloatType`: it is process-wide and silently changes the
  precision of every model built afterwards. `ECOSMO` takes `FT` from the grid
  (`EcosmoParameters(FT)`).
- Never extend `getproperty` to fix an undefined-property bug — fix the caller.
- A variable named the same as a function produces "type is not callable" — rename the variable.
- Never add, remove or change `[deps]` in `Project.toml` (or `box_model/Project.toml`) on your own
  initiative; touch `[compat]` only when asked. If a dependency would meaningfully simplify the
  implementation, ask — including packages already present transitively.
- `Pkg.resolve()` only checks the manifest against `[compat]`; after raising a bound use
  `Pkg.update()`.
- `Base.@kwdef` on a parametric struct does not convert its arguments:
  `QuadraticBottomStress(reference_density = 1026)` is a `MethodError` (an `Int` next to the
  `Float64` default); write `1026.0`. `ecosmo_sediment` converts it to the grid's `FT` afterwards.
