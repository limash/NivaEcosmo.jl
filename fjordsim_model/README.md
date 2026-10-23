# FjordSim setups with ECOSMO

[FjordSim](https://github.com/NIVANorge/FjordSim.jl) setups that run the NivaEcosmo biogeochemistry.

| File | Content |
|---|---|
| `moro.jl` | Møre og Romsdal, 5.12-9.00°E by 62.00-63.43°N at ~500 m, 59 levels to 1137 m, the year 2023 |
| `Project.toml` | the environment: FjordSim from GitHub, NivaEcosmo from this repository |
| `../data/fabm.yaml` | the ECOSMO parameters (A20 v3aa), the carbonate pH scale and the sediment's initial pools |
| `../data/moro_ecosmo.yaml` | Moro's boundary, initial and river values of the ECOSMO tracers, with their seasonal cycle. Placeholders |

## Install

Julia 1.12 or later, since FjordSim is run with `julia -m`.

```bash
cd NivaEcosmo.jl
julia --project=fjordsim_model -e 'using Pkg; Pkg.instantiate()'
```

`Pkg.instantiate` clones FjordSim's `main` branch from GitHub and uses NivaEcosmo from `..`, so an
edit to NivaEcosmo takes effect in the next run. To run against a local FjordSim checkout instead:

```bash
julia --project=fjordsim_model -e 'using Pkg; Pkg.develop(path = "path/to/FjordSim.jl")'
```

Check that the setup loads:

```bash
julia --project=fjordsim_model -e 'using FjordSim; fjord_config("fjordsim_model/moro.jl")'
```

Before running the steps:

- `add_rivers` reads NVE's HydAPI key from `NVE_API_KEY` (free at https://hydapi.nve.no/Users).
- Data goes to `~/FjordSim_data/moro/` and results to `~/FjordSim_results/moro/`, unless
  `FJORDSIM_DATA_ROOT` and `FJORDSIM_RESULTS_ROOT` point elsewhere.
- `prepare_bathymetry` downloads the Geonorge FileGDB (4.5 GB) into `~/FjordSim_data/moro/geonorge/`
  if it is not there. A copy from another FjordSim setup can be put there instead.
- `run_simulation` needs a CUDA GPU.

## Run

Each step is

```bash
julia --project=fjordsim_model -m FjordSim <step> --config fjordsim_model/moro.jl
```

in this order: `prepare_bathymetry`, `download_initial_state`, `prepare_initial_state`,
`add_rivers`, `download_boundaries`, `prepare_boundaries`, `download_atmosphere`,
`prepare_atmosphere`, `run_simulation`, `validate_simulation`.

A changed `data/moro_ecosmo.yaml` needs `add_rivers` and `prepare_boundaries` re-run. FjordSim's
`docs/running.md` lists what else to re-run after a change.

FjordSim's `gcp/fjordsim-gcp` runs a setup in FjordSim's own environment, which has no NivaEcosmo,
so these setups run locally for now.

## What to check while preparing the data

Each step writes its file, and most a plot beside it, to `~/FjordSim_data/moro/`. Look at them
before going on: a problem found later means re-running every step after the one that caused it.

**After `prepare_bathymetry`** (`bathymetry.png`, `bathymetry.nc`):

- **Deepest cell against the deepest level.** The deepest cell should be close to the first entry
  of `z_faces`. Deeper water than that is cut off at the deepest level. That is harmless for a few
  columns near an open edge, but otherwise add a level. A deepest level with no water in it wastes
  time, so remove it.
- **Open edges.** An edge with water along it must be in `open_edges`, or it is a wall. An edge
  that is land, or touches water in only a cell or two, stays closed.
- **Invented water.** Water where there is land, inland or up a valley, means the depth data had a
  gap that the gridding filled. Set `max_sounding_distance` (2000 m worked for Oslofjord) and re-run.
- **Sounds and sills.** Narrow sounds between islands should still be open, and known sills close
  to their real depth.
- **Pits on open edges.** An edge cell much deeper than its neighbours traps water. On the Oslofjord
  run such cells sat at the freezing point. Note them and watch them in the first snapshots.

These numbers are quick to read from the file:

```julia
using NCDatasets
h = NCDataset(ds -> ds["h"][:, :], expanduser("~/FjordSim_data/moro/bathymetry.nc"))
minimum(h)                                            # deepest cell, m (negative)
count(<(0), h[1, :]) / size(h, 2)                     # western edge: fraction that is water
count(<(0), h[:, end]) / size(h, 1)                   # northern edge
```

For the current Moro bathymetry: the deepest cell is 1153 m against a deepest level at 1137 m,
which cuts 30 slope columns inside the sponge by up to 16 m; the western and northern edges are
85-86 % water and open, the southern and eastern ones land but for one fjord-head cell each and
closed. There is no invented water. The edge pits are on the western edge at 62.25°N and the
northern edge at 8.71°E and 8.73°E.

**After `add_rivers`** (its log, `rivers_nve.png`): the log lists every river mouth it found, with
its `vassdragsnr`, discharge and the number of columns it enters through. Check that the large
rivers (Driva, Rauma, Surna, Eira, Aura) sit at the right fjord heads. Without a gauge a river runs
at its catchment's mean flow all year, so attach their HydAPI gauges as `NVERiver` entries in
`outlets`, keyed by that `vassdragsnr` (FjordSim's `examples/lista.jl` shows how), and re-run
`add_rivers`.

**After `prepare_boundaries`** (`boundaries.png`): the open edges carry data wherever there is
water, and the log shows the ECOSMO tracers being appended (`Writing west_no3 ...`).

**After `prepare_atmosphere`** (`atmosphere.png`): the fields cover the whole box and the whole
year.
