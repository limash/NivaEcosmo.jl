# ECOSMO as an OceanBioME biogeochemistry: the tracer tendencies of `interior_rates`, the sinking
# of `sinking_speeds`, and the diagnostics the tendencies depend on — day length and the carbonate
# system — kept as auxiliary fields and refreshed in `update_biogeochemical_state!`, which is when
# FABM computes them too: from the current state, before any rate is evaluated.

# The water-column state variables, in FABM's order. `dic` and `alk` belong to the carbonate
# module in the Fortran; here they are ECOSMO's own, which is what `couple_co2: true` amounts to.
const ECOSMO_TRACERS = (:no3, :nh4, :pho, :sil, :oxy, :fla, :dia, :flachl, :diachl,
                        :microzoo, :mesozoo, :det, :opa, :cal, :dom, :dic, :alk)

# The ones whose `interior_rates` input is the concentration itself (`dic` and `alk` act through
# the carbonate diagnostics instead).
const ECOSMO_REACTIVE_TRACERS = ECOSMO_TRACERS[1:15]

const ECOSMO_SINKING_TRACERS = (:det, :opa, :cal)

# FABM's `minimum` for each state variable, which `check_state(repair = .true.)` clips to.
const ECOSMO_MINIMA = (no3 = 0.0, nh4 = 0.0, pho = 0.0, sil = 0.0, oxy = 0.0,
                       fla = 1.0e-7, dia = 1.0e-7, flachl = 1.0e-7 / 20.0, diachl = 1.0e-7 / 27.0,
                       microzoo = 1.0e-7, mesozoo = 1.0e-7, det = 0.0, opa = 0.0, cal = 0.0, dom = 0.0,
                       dic = 0.0, alk = 1.0e-4)

# The carbonate diagnostics kept between steps. Every one of them is needed, not just reported:
# when the carbonate iteration fails the Fortran keeps the previous pH and pCO2 and computes the
# saturation state from the previous carbonate ion concentration.
const CARBONATE_FIELDS = (:pH, :pCO2, :Carb, :Om_cal)

"""
    EcosmoBiogeochemistry

The OceanBioME underlying biogeochemistry for ECOSMO. Built by [`ECOSMO`](@ref), which wraps it
in an OceanBioME `Biogeochemistry` together with its light and sediment.
"""
struct EcosmoBiogeochemistry{FT, K, W, C, D, L, R} <: AbstractBiogeochemistry
    parameters :: EcosmoParameters{FT}
    constants :: K
    pH_scale :: Int
    sinking_velocities :: W
    carbonate :: C
    day_length :: D
    latitude :: L
    start_day_of_year :: FT
    density :: R
end

Adapt.adapt_structure(to, bgc::EcosmoBiogeochemistry) =
    EcosmoBiogeochemistry(bgc.parameters,
                          bgc.constants,
                          bgc.pH_scale,
                          adapt(to, bgc.sinking_velocities),
                          adapt(to, bgc.carbonate),
                          adapt(to, bgc.day_length),
                          bgc.latitude,
                          bgc.start_day_of_year,
                          bgc.density)

"""
    ECOSMO(grid;
           parameters = EcosmoParameters(eltype(grid)),
           surface_shortwave = 0,
           light_attenuation = GOTMLight(grid, surface_shortwave),
           sediment = ecosmo_sediment(grid; parameters),
           latitude = nothing,
           start_day_of_year = 0,
           pH_scale = 1,
           density = teos10_polynomial_approximation,
           repair_state = true,
           modifiers = nothing)

The NIVA ECOSMO model (`niva/niva_ecosmo` with `niva_ecosmo_daylength`, `niva_ecosmo_carbonate`
and `gotm/light`) as an OceanBioME `Biogeochemistry`, to pass to an Oceananigans model as
`biogeochemistry`.

The tracers keep their FABM names and units (`ECOSMO_TRACERS`): nutrients, plankton and organic
matter in mgC m⁻³ (Redfield), chlorophyll in mg m⁻³, oxygen, `dic` and `alk` in mmol m⁻³. They are
not required to include `T` and `S`, but the model must have both, in °C and PSU.

Keyword arguments
=================

- `parameters`: the [`EcosmoParameters`](@ref); [`fabm_yaml_parameters`](@ref) reads them from a
  `fabm.yaml`.
- `surface_shortwave`: downwelling shortwave at the surface (W m⁻²) for the default light model.
- `light_attenuation`: the light model; it must provide the `PAR` auxiliary field in W m⁻².
- `sediment`: the benthic model, an [`ecosmo_sediment`](@ref) by default. `nothing` leaves the
  bottom closed, so sinking material accumulates in the bottom cell.
- `latitude`: degrees north for the day length; `nothing` reads it from a latitude-longitude grid.
- `start_day_of_year`: the day of the year at model time zero, counted from 0 at midnight on
  1 January (so 1 January 06:00 is `0.25`).
- `pH_scale`: the carbonate module's `pHscale`: 1 total, 0 seawater, -1 old seawater formulation.
- `density`: `density(T, S, P)` in kg m⁻³ with `P` in bar, for the carbonate unit conversions.
- `repair_state`: clip every tracer to its FABM minimum before each tendency evaluation, as a FABM
  host does with `repair_state`. The equations divide by phytoplankton biomass and need it.
- `modifiers`: further OceanBioME modifiers, applied after the repair.
"""
function ECOSMO(grid::AbstractGrid{FT};
                parameters = EcosmoParameters(FT),
                surface_shortwave = 0,
                light_attenuation = GOTMLight(grid, surface_shortwave),
                sediment = ecosmo_sediment(grid; parameters),
                latitude = nothing,
                start_day_of_year = 0,
                pH_scale = 1,
                density = teos10_polynomial_approximation,
                repair_state = true,
                modifiers = nothing) where FT

    pH_scale in (-1, 0, 1) || throw(ArgumentError("pH_scale must be -1, 0 or 1, got $pH_scale"))
    isnothing(latitude) && !(grid isa LatitudeLongitudeGrid ||
                             (grid isa ImmersedBoundaryGrid && grid.underlying_grid isa LatitudeLongitudeGrid)) &&
        throw(ArgumentError("`latitude` must be given for a grid without latitude"))

    carbonate = NamedTuple{CARBONATE_FIELDS}(map(_ -> CenterField(grid), CARBONATE_FIELDS))
    # FABM's `missing_value` for Om_cal; the others are missing as zero, which `CenterField` is.
    set!(carbonate.Om_cal, 4)

    underlying = EcosmoBiogeochemistry(parameters,
                                       ecosmo_constants(FT),
                                       Int(pH_scale),
                                       sinking_velocity_fields(grid, parameters),
                                       carbonate,
                                       Field{Center, Center, Nothing}(grid),
                                       isnothing(latitude) ? nothing : convert(FT, latitude),
                                       convert(FT, start_day_of_year),
                                       density)

    if repair_state
        repair = ClipToMinimum(map(value -> convert(FT, value), ECOSMO_MINIMA))
        modifiers = isnothing(modifiers) ? repair :
                    modifiers isa Tuple ? (repair, modifiers...) : (repair, modifiers)
    end

    return Biogeochemistry(underlying; light_attenuation, sediment, modifiers)
end

"""
    sinking_velocity_fields(grid, parameters)

The vertical drift of detritus, opal and calcite on `grid`'s z-faces, downward at
`sinking_speeds`' speed for the face's depth.

Zero on the seabed and at the surface: as for any FABM host, nothing sinks out of the water column
by advection. What reaches the bottom cell leaves it through the sediment's settling flux
(`bottom_fluxes`).

The seabed includes an immersed one. Oceananigans does not close an immersed face to a
biogeochemical drift velocity, and an immersed cell's tracers are zeroed on every state update, so
without this whatever sank through it would be lost.
"""
function sinking_velocity_fields(grid, parameters)
    return NamedTuple{ECOSMO_SINKING_TRACERS}(map(ECOSMO_SINKING_TRACERS) do name
        w = ZFaceField(grid)
        launch!(architecture(grid), grid, :xyz, _set_sinking_velocity!, w, grid, parameters, Val(name))
        w
    end)
end

@kernel function _set_sinking_velocity!(w, grid, parameters, ::Val{name}) where name
    i, j, k = @index(Global, NTuple)

    depth = -znode(i, j, k, grid, Center(), Center(), Face())
    speed = getproperty(sinking_speeds(parameters, depth), name)
    seabed = (k == 1) | immersed_cell(i, j, k - 1, grid)

    @inbounds w[i, j, k] = ifelse(seabed, zero(speed), -speed)
end

required_biogeochemical_tracers(::EcosmoBiogeochemistry) = ECOSMO_TRACERS

biogeochemical_auxiliary_fields(bgc::EcosmoBiogeochemistry) =
    merge(bgc.carbonate, (; day_length = bgc.day_length))

for name in ECOSMO_SINKING_TRACERS
    @eval biogeochemical_drift_velocity(bgc::EcosmoBiogeochemistry, ::Val{$(QuoteNode(name))}) =
        (u = ZeroField(), v = ZeroField(), w = bgc.sinking_velocities.$name)
end

@inline function cell_rates(i, j, k, bgc::EcosmoBiogeochemistry, fields, auxiliary_fields)
    c = @inbounds NamedTuple{ECOSMO_REACTIVE_TRACERS}(map(name -> fields[name][i, j, k], ECOSMO_REACTIVE_TRACERS))

    temp = @inbounds fields.T[i, j, k]
    D = @inbounds auxiliary_fields.day_length[i, j, 1]
    E0 = @inbounds auxiliary_fields.PAR[i, j, k]
    om_cal = @inbounds auxiliary_fields.Om_cal[i, j, k]

    return interior_rates(bgc.parameters, bgc.constants, c, temp, D, E0, om_cal)
end

# Each tracer evaluates the whole model and keeps its own rate; what it does not use is compiled
# away. Anything else the model carries (`T`, `S`) has no biogeochemical source.
@inline (bgc::EcosmoBiogeochemistry)(i, j, k, grid, val_name, clock, fields, auxiliary_fields) = zero(grid)

for name in ECOSMO_TRACERS
    @eval @inline (bgc::EcosmoBiogeochemistry)(i, j, k, grid, ::Val{$(QuoteNode(name))}, clock, fields, auxiliary_fields) =
        cell_rates(i, j, k, bgc, fields, auxiliary_fields).rates.$name
end

"""
    day_of_year(bgc, clock)

FABM's `number_of_days_since_start_of_the_year` at `clock`, on a 365-day calendar.
"""
@inline day_of_year(bgc::EcosmoBiogeochemistry, clock) = mod(bgc.start_day_of_year + clock.time / SECONDS_PER_DAY, 365)

@inline column_latitude(i, j, grid, latitude) = latitude
@inline column_latitude(i, j, grid, ::Nothing) = φnode(i, j, 1, grid, Center(), Center(), Center())

function update_biogeochemical_state!(model, bgc::EcosmoBiogeochemistry)
    grid = model.grid
    arch = architecture(grid)

    launch!(arch, grid, :xy, _compute_day_length!,
            bgc.day_length, grid, bgc.constants, bgc.latitude, day_of_year(bgc, model.clock))

    tracers = model.tracers
    launch!(arch, grid, :xyz, _compute_carbonate!,
            bgc.carbonate, grid, bgc.constants, tracers.T, tracers.S, tracers.dic, tracers.alk, bgc.density,
            bgc.pH_scale)

    return nothing
end

@kernel function _compute_day_length!(D, grid, constants, latitude, yday)
    i, j = @index(Global, NTuple)
    @inbounds D[i, j, 1] = day_length(constants, column_latitude(i, j, grid, latitude), yday)
end

@kernel function _compute_carbonate!(carbonate, grid, constants, T, S, dic, alk, density, pH_scale)
    i, j, k = @index(Global, NTuple)

    @inbounds if !immersed_cell(i, j, k, grid)
        pressure = -znode(i, j, k, grid, Center(), Center(), Center())  # dbar ≈ m
        ρ = density(T[i, j, k], S[i, j, k], pressure / 10)

        previous = (pH = carbonate.pH[i, j, k], pCO2 = carbonate.pCO2[i, j, k],
                    CarbA = zero(ρ), BiCarb = zero(ρ), Carb = carbonate.Carb[i, j, k])

        d = carbonate_diagnostics(constants, dic[i, j, k], alk[i, j, k], T[i, j, k], S[i, j, k], ρ, pressure,
                                  previous, pH_scale)

        carbonate.pH[i, j, k] = d.pH
        carbonate.pCO2[i, j, k] = d.pCO2
        carbonate.Carb[i, j, k] = d.Carb
        carbonate.Om_cal[i, j, k] = d.Om_cal
    end
end

Base.summary(::EcosmoBiogeochemistry{FT}) where {FT} = "NIVA ECOSMO{$FT}"
Base.show(io::IO, bgc::EcosmoBiogeochemistry) =
    print(io, summary(bgc), " with tracers ", ECOSMO_TRACERS)

#####
##### FABM's state repair
#####

"""
    ClipToMinimum(minima::NamedTuple)

An OceanBioME modifier that clips each tracer named in `minima` to its minimum in place before the
tendencies are computed — what a FABM host's `check_state` does with `repair = .true.`.

Clipping is not mass conserving. It is here because FABM's ECOSMO relies on it: the chlorophyll
equations divide by phytoplankton carbon, which must stay above zero.
"""
struct ClipToMinimum{M}
    minima :: M
end

function update_biogeochemical_state!(model, clip::ClipToMinimum)
    for (name, minimum) in pairs(clip.minima)
        tracer = parent(model.tracers[name])
        tracer .= max.(minimum, tracer)
    end
    return nothing
end

Base.summary(clip::ClipToMinimum) = "ClipToMinimum on $(keys(clip.minima))"
