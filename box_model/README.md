# Drammensfjorden box model (Julia ECOSMO)

The Julia twin of niva-ecosmo's `box_model` (`run_box.py`, which runs the Fortran model through
pyfabm). A single well-mixed box with prescribed temperature, salinity, light and wind, a sediment
layer underneath, and exchange with the river and the outer fjord.

It exists to show that the Julia port computes what the Fortran computes. The driver repeats
`run_box.py` step for step: forward Euler with `DT = 600 s`, clipping every state variable to its
FABM minimum after each step (`check_state(repair=True)`), surface and bottom fluxes divided by the
box depth, the Knudsen exchange, the diel light cycle, and numpy's interpolation of the daily
forcing.

## Files

| File | Content |
|---|---|
| `../data/fabm.yaml`, `../data/forcing.csv`, `../data/boundary.yaml` | the inputs, copied unchanged from niva-ecosmo's `box_model/` |
| `run_box.jl` | the box model; writes `output.nc` and `output.png` |
| `compare.jl` | compares `output.nc` with the pyfabm `output.nc`, writes `compare.png` |
| `Project.toml` | the environment; uses the package from `..` |

`fabm.yaml` is read with `NivaEcosmo.fabm_configuration`, which finds the `niva/niva_ecosmo`,
`niva/niva_ecosmo_carbonate` and `gotm/light` instances and takes their parameters and
initialization from the file, exactly as FABM does.

## Run

```bash
cd ~/src/NivaEcosmo.jl
julia --project=box_model -e 'using Pkg; Pkg.instantiate()'    # once
julia --project=box_model box_model/run_box.jl                 # full run -> output.nc, output.png
julia --project=box_model box_model/run_box.jl --closed        # no exchange, no burial: total P is conserved
julia --project=box_model box_model/compare.jl ~/src/niva-ecosmo/box_model/output.nc
```

`run_box.jl` takes `--yaml`, `--forcing`, `--boundary` and `--out` like `run_box.py`, and
`--no-plot`.

## Agreement with the Fortran

Against `~/src/niva-ecosmo/box_model/output.nc` (pyfabm, niva-ecosmo commit `b351e6f`), the largest
difference over the year in any state variable or diagnostic is 3.3 × 10⁻¹⁵ of that variable's
largest value, and the `--closed` run conserves total phosphorus to the same relative 1.8 × 10⁻¹⁴
as pyfabm's. Step by step, given the same state, the rates agree to a few units in the last place
and most are bit-identical.

What remains is the rounding of `exp`, `log`, `pow` and the trigonometric functions, which Julia
implements itself and gfortran takes from glibc. Everything else rounds the same: the expressions
keep the Fortran's evaluation order, `x**37` is the multiplication chain gfortran emits, and the
two single-precision literals of the carbonate module are single precision here too.
