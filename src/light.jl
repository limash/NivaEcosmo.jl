# FABM's `gotm/light`, the light model the ECOSMO configurations run with: shortwave split into a
# non-visible fraction `A`, absorbed over `g1`, and PAR, attenuated by the water over `g2` and by
# chlorophyll. ECOSMO reads its `E0` from the PAR at the layer centre.

"""
    chlorophyll_attenuation_weight(Exphy, ExGroup)

The weight a phytoplankton group's chlorophyll carries in the PAR attenuation coefficient, as FABM
actually assembles it.

`niva_ecosmo` registers each chlorophyll variable with the attenuation aggregate twice, once with
`Exphy` and once with the group's own `ExPs` or `ExPl`, intending the two to add. FABM keeps one
contribution per variable and aggregate, though: a zero factor is not registered at all, and a
second non-zero one *replaces* the first (`contribution_list_add` in `fabm_types.F90`). The group
coefficient therefore wins whenever it is non-zero. This reproduces that, not the comment's intent.
"""
@inline chlorophyll_attenuation_weight(Exphy, ExGroup) = ifelse(ExGroup != 0, ExGroup, Exphy)

"""
    attenuation_coefficient(p::EcosmoParameters, flachl, diachl)

The chlorophyll contribution to the PAR attenuation coefficient, m⁻¹.
"""
@inline attenuation_coefficient(p::EcosmoParameters, flachl, diachl) =
    chlorophyll_attenuation_weight(p.Exphy, p.ExPs) * flachl +
    chlorophyll_attenuation_weight(p.Exphy, p.ExPl) * diachl

"""
    gotm_light(swr0, A, g1, g2, z, bioext)

PAR and shortwave (W m⁻²) at depth `z` below a surface shortwave flux `swr0`, where `bioext` is
the chlorophyll attenuation integrated from the surface to `z`: `(; par, swr)`.
"""
@inline function gotm_light(swr0, A, g1, g2, z, bioext)
    par = swr0 * (1.0 - A) * exp(-z / g2 - bioext)
    swr = par + swr0 * A * exp(-z / g1)
    return (; par, swr)
end

"""
    GOTMLight(grid, surface_shortwave; A = 0.58, g1 = 0.35, g2 = 23.0,
              parameters = nothing, discrete_form = false)

An OceanBioME light attenuation component that computes PAR as FABM's `gotm/light` does, with the
chlorophyll attenuation of the ECOSMO model it is part of. Exposes the cell-centred PAR (W m⁻²) as
the auxiliary field `PAR`, which ECOSMO reads as its `E0`.

`surface_shortwave` is the downwelling shortwave just below the surface, in W m⁻²: a number, a
field, or a function `f(x, y, t)` (with `parameters`, `f(x, y, t, parameters)`; with
`discrete_form = true`, `f(i, j, grid, clock, fields)`), as for OceanBioME's own light models.

`A` is the non-visible fraction of shortwave, `g1` its e-folding depth (m) and `g2` the e-folding
depth of PAR in water without chlorophyll (m) — the `gotm/light` parameters `a`, `g1` and `g2`.
"""
struct GOTMLight{FT, F, S}
    A :: FT
    g1 :: FT
    g2 :: FT
    field :: F
    surface_shortwave :: S
end

function GOTMLight(grid::AbstractGrid{FT}, surface_shortwave;
                   A = 0.58, g1 = 0.35, g2 = 23.0,
                   parameters = nothing, discrete_form = false) where FT
    # OceanBioME's helper for its PAR field also regularizes the surface condition into something
    # `getbc` can evaluate. It gives the field's top halo the shortwave value, which nothing reads.
    field, surface_shortwave = PAR_field(grid, surface_shortwave, parameters, discrete_form)
    return GOTMLight(convert(FT, A), convert(FT, g1), convert(FT, g2), field, surface_shortwave)
end

Adapt.adapt_structure(to, light::GOTMLight) =
    GOTMLight(light.A, light.g1, light.g2, adapt(to, light.field), adapt(to, light.surface_shortwave))

biogeochemical_auxiliary_fields(light::GOTMLight) = (PAR = light.field,)

function update_biogeochemical_state!(model, light::GOTMLight)
    p = model.biogeochemistry.underlying_biogeochemistry.parameters
    grid = model.grid

    launch!(architecture(grid), grid, :xy, _integrate_gotm_light!,
            light.field, grid, model.clock, fields(model), light.surface_shortwave,
            light.A, light.g1, light.g2, p, model.tracers.flachl, model.tracers.diachl)

    return nothing
end

@kernel function _integrate_gotm_light!(PAR, grid, clock, model_fields, surface_shortwave,
                                        A, g1, g2, p, flachl, diachl)
    i, j = @index(Global, NTuple)

    swr0 = getbc(surface_shortwave, i, j, grid, clock, model_fields)
    z = zero(grid)
    bioext = zero(grid)

    @inbounds for k in grid.Nz:-1:1
        dz = Δzᶜᶜᶜ(i, j, k, grid)
        ext = attenuation_coefficient(p, flachl[i, j, k], diachl[i, j, k])

        z = z + dz * 0.5
        bioext = bioext + ext * dz * 0.5

        PAR[i, j, k] = gotm_light(swr0, A, g1, g2, z, bioext).par

        z = z + dz * 0.5
        bioext = bioext + ext * dz * 0.5
    end
end

Base.summary(::GOTMLight{FT}) where {FT} = "GOTMLight{$FT}"
Base.show(io::IO, light::GOTMLight) =
    print(io, summary(light), " (A = $(light.A), g1 = $(light.g1) m, g2 = $(light.g2) m)")
