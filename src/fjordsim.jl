# The FjordSim components of ECOSMO, implemented in `ext/NivaEcosmoFjordSimExt.jl`, which Julia loads
# when a setup loads FjordSim and NCDatasets beside NivaEcosmo. A tracer value is
# `(z, date, source) -> value`: `z` the depth (m, negative down), `date` a `DateTime`, and `source` a
# `NamedTuple` of the source dataset's prepared values at that point, empty without one.

"""
    ecosmo_biogeochemistry_config(; configuration, surface_shortwave, bottom_drag_coefficient, start_day_of_year)

A FjordSim biogeochemistry config building `ECOSMO` on the simulation grid, with its sediment pools at
`configuration`'s (`fabm_configuration`'s) initial values. `bottom_drag_coefficient` is the sediment's
quadratic bottom stress, which should be the ocean's own drag.
"""
function ecosmo_biogeochemistry_config end

"""
    ecosmo_gas_exchange(; wind_speed, atmospheric_pCO2, pH_scale)

A FjordSim boundary-condition config giving `oxy` and `dic` ECOSMO's air-sea O₂ and CO₂ exchange as
top conditions. `pH_scale` must be the one `ECOSMO` gets.
"""
function ecosmo_gas_exchange end

"""
    ecosmo_boundaries(base; sources = nothing, tracers)

A FjordSim open-boundary config: `base`'s prepared boundary file, with `<edge>_<tracer>` appended for
every entry of `tracers`, a `NamedTuple` of `(z, date, source) -> value`. `sources`, a second boundary
config over the same edges and period (TOPAZ5-ECOSMO's), is downloaded and prepared beside `base` and
gives each point's `source`.
"""
function ecosmo_boundaries end

"""
    ecosmo_initial_state(base; sources = nothing, tracers)

A FjordSim initial-state config: `base`'s prepared state with every entry of `tracers` appended,
evaluated in each water cell at `base`'s `start_date`, so `FromInitialState` starts the tracers from
it. `sources`, a second initial-state config, gives each cell's `source`.
"""
function ecosmo_initial_state end

"""
    ecosmo_rivers(base, concentrations)

A FjordSim river config: `base`'s rivers, each also carrying every entry of `concentrations`, a
`NamedTuple` of `date -> value`.
"""
function ecosmo_rivers end

"""
    ecosmo_surface_writer(; name, output_file, interval, split_times, overwrite_existing)

A FjordSim writer of the top cell of T, S, u, v, e, every ECOSMO tracer, the carbonate system and
`PAR`, of `η` and `day_length`, the sediment pools, and the O₂ and CO₂ air-sea fluxes (positive into
the water), every `interval` seconds, starting a new file part at each of `split_times`.
"""
function ecosmo_surface_writer end
