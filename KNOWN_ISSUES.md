# Known issues

Found while porting niva-ecosmo (commit `b351e6f`). NivaEcosmo copies the Fortran's behaviour,
including issues 1–3, so that the two implementations give the same numbers. Issues 4 and 5 are in
the Julia packages, and NivaEcosmo works around them.

| # | Issue | Where | Affects the shipped `fabm.yaml`? |
|---|---|---|---|
| 1 | `Exphy` and `ExPs`/`ExPl` replace each other instead of adding (**needs a decision**) | Fortran + FABM | no |
| 2 | `TSOM_s` counts only phosphorus | Fortran + FABM | no (not used) |
| 3 | Settling uses the minimum sinking speed | Fortran | **yes, in 3D** |
| 4 | Sinking through an immersed seabed loses mass | Oceananigans | fixed in NivaEcosmo |
| 5 | Sediment fluxes are divided by the wrong cell height | OceanBioME | fixed in NivaEcosmo |

The background for issues 1 and 2: an *aggregate variable* in FABM is a sum that the models add
to, such as the total light attenuation. FABM stores one weight per variable in each sum. If a
model adds the same variable a second time, the new weight replaces the old one, unless it is zero:
a zero weight is skipped (`contribution_list_add` in FABM's `fabm_types.F90`).

## 1. `Exphy` and `ExPs`/`ExPl` replace each other (decision needed)

The Fortran adds the chlorophyll of small phytoplankton (`flachl`) to the light attenuation twice:
once with `Exphy`, the coefficient for all chlorophyll, and once with `ExPs`, the coefficient for
small phytoplankton. `diachl` is treated the same way with `ExPl` (`niva_ecosmo.F90`, lines
638–646). The comment of 09/06/2026 says the two should add up. Because of the rule above, FABM
uses `ExPs` alone whenever `ExPs` is not zero.

**Example.** `Exphy = 0.041` (as in `fabm.yaml`), `ExPs = 0.02`, `flachl = 2 mg m⁻³`:

- intended: (0.041 + 0.02) × 2 = 0.122 m⁻¹
- FABM and NivaEcosmo: 0.02 × 2 = 0.04 m⁻¹, so the water is three times clearer than intended.

**Today.** `fabm.yaml` does not set `ExPs`/`ExPl`. They default to 0, so `Exphy` is used, as
intended.

**Decision.** Keep FABM's behaviour, or make the two add? To make them add:

- in Julia, change `chlorophyll_attenuation_weight` in `src/light.jl` to `Exphy + ExGroup`;
- in the Fortran, add `flachl` once with the weight `Exphy + ExPs` (and `diachl` with
  `Exphy + ExPl`).

Change both, or the two implementations stop matching.

## 2. `TSOM_s` counts only phosphorus

`TSOM_s` (total suspended organic mass, mg L⁻¹) is meant to add up the carbon, nitrogen and
phosphorus mass of each organic pool. The Fortran adds each pool three times, with the C, N and P
weights (lines 654 onwards). By the rule above, only the last weight, phosphorus, is kept.

**Example.** `fla = 100 mgC m⁻³`:

- intended: 100 × 10⁻³ × 1.5547 × (1 + 0.176 + 0.024) = 0.187 mg L⁻¹
- FABM: 100 × 10⁻³ × 1.5547 × 0.024 = 0.0038 mg L⁻¹, about 50 times too small.

**Today.** Only `light_iop_aksnes.F90` reads `TSOM_s`, and it is not part of the niva-ecosmo build,
so nothing is affected. NivaEcosmo does not port `TSOM_s`/`TSIM_s`. Fix it in the Fortran by adding
each pool once with the summed weight.

## 3. Settling uses the minimum sinking speed

Detritus sinks faster with depth: `sinkD + cz_sinkD × depth` (`get_vertical_movement`, line 1263).
It settles from the bottom cell into the sediment at `sinkD`, the speed at the surface
(`do_bottom`, lines 1176–1182). Opal and calcite are treated the same way. Note 9a in the Fortran
sets the settling rate equal to the sinking speed. It was written before the depth dependence was
added on 21/05/2026.

**Example.** `fabm.yaml` has `sinkD = 5 m d⁻¹` and `cz_sinkD = 0.039 d⁻¹`. At a 200 m seabed,
detritus arrives in the bottom cell at about 5 + 0.039 × 200 ≈ 13 m d⁻¹ but leaves it at
5 m d⁻¹. It piles up until the bottom cell holds about 13 / 5 ≈ 2.6 times the concentration of
the cell above. More detritus is then remineralized just above the seabed, and less reaches the
sediment.

**Today.** This affects 3D runs with the shipped `fabm.yaml`. It does not affect the box model,
which has no water column to sink through. If this is not intended, settle at the sinking speed at
the seabed: in `do_bottom` in the Fortran, and in `bottom_fluxes` (`src/ecosmo.jl`) in Julia.

## 4. Oceananigans: sinking through an immersed seabed loses mass

An *immersed* seabed is how Oceananigans represents bathymetry: the cells below the seabed are
"land" inside the grid. Oceananigans (0.113) lets particles sink through the seabed into those
land cells, then sets every land cell's tracers to zero at each step. Whatever sinks out of the
deepest water cell is therefore deleted.

**Example.** In a 4 × 3 × 6 test grid with a varying seabed, 20 steps lost 0.1–0.3 % of the
phosphorus.

**Fix in NivaEcosmo.** The sinking speed is zero on every face above a seabed cell, immersed or
not, as in FABM. Material leaves the bottom cell only through the sediment's settling flux
(`sinking_velocity_fields` in `src/biogeochemistry.jl`). The same test now conserves phosphorus to
about 10⁻¹⁶ (`test/runtests.jl`, "3D run conserves phosphorus").

## 5. OceanBioME: sediment fluxes are divided by the wrong cell height

A flux between sediment and water (mmol m⁻² d⁻¹) becomes a concentration change in the bottom cell
(mmol m⁻³ d⁻¹) when it is divided by that cell's height. OceanBioME (0.20) divides by the distance
between the bottom cell's centre and the centre of the cell below (`Δzᶜᶜᶠ`,
`Sediments/tracer_coupling.jl`, line 37). It should divide by the bottom cell's height (`Δzᶜᶜᶜ`).
The two differ whenever the cells around the seabed have different heights: with an immersed
seabed on a stretched grid, or with partial cells.

**Example.** The seabed is immersed at 10 m depth. The bottom water cell is 6 m high, and the land
cell below it is 2 m high. Their centres are 4 m apart. The sediment releases 12 mmol m⁻² d⁻¹:

- correct: 12 / 6 = 2 mmol m⁻³ d⁻¹, so the water column gains 12 mmol m⁻² d⁻¹
- OceanBioME: 12 / 4 = 3 mmol m⁻³ d⁻¹, so the water column gains 18 mmol m⁻² d⁻¹, 50 % more than
  the sediment lost.

**Fix in NivaEcosmo.** `src/sediment.jl` replaces OceanBioME's `update_tendencies!` for the ECOSMO
sediment with one that divides by `Δzᶜᶜᶜ`. Worth reporting to OceanBioME.
