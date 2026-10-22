# ECOSMO's benthic layer (`do_bottom` and the bottom state variables `sed1`–`sed4`) as an OceanBioME
# sediment biogeochemistry. OceanBioME's `BiogeochemicalSediment` holds the four pools as flat
# fields, steps them with the ocean model's time stepper, copies the bottom-cell values the fluxes
# need into `tracked_fields`, and adds the fluxes to the bottom cell's tendencies.

const ECOSMO_SEDIMENT_FIELDS = (:sed1, :sed2, :sed3, :sed4)

# The tracers `do_bottom` exchanges with the water.
const ECOSMO_BOTTOM_COUPLED_TRACERS = (:det, :no3, :nh4, :opa, :sil, :pho, :cal, :oxy, :dic, :alk)

# What `do_bottom` reads from the bottom cell: tracers, temperature and the calcite saturation.
const ECOSMO_BOTTOM_TRACKED_FIELDS = (:T, :oxy, :no3, :det, :opa, :cal, :Om_cal)

"""
    QuadraticBottomStress(; drag_coefficient = 0.003, reference_density = 1026.0)

Bottom stress (Pa) for ECOSMO's resuspension threshold from the bottom-cell velocity,
ρ₀ C_d (u² + v²). In FABM the host supplies `bottom_stress`; this is the host's job done here, and
the coefficient should be the one the ocean model's own bottom drag uses.
"""
Base.@kwdef struct QuadraticBottomStress{FT}
    drag_coefficient :: FT = 0.003
    reference_density :: FT = 1026.0
end

"""
    EcosmoSediment(parameters, bottom_stress)

The sediment biogeochemistry of ECOSMO, for `OceanBioME.Sediments.BiogeochemicalSediment`. Build
the coupled sediment with [`ecosmo_sediment`](@ref).

`bottom_stress` is a constant in Pa or a [`QuadraticBottomStress`](@ref).
"""
struct EcosmoSediment{FT, B} <: AbstractSedimentBiogeochemistry
    parameters :: EcosmoParameters{FT}
    bottom_stress :: B
end

"""
    ecosmo_sediment(grid; parameters = EcosmoParameters(eltype(grid)), bottom_stress = 0)

ECOSMO's sediment on `grid`, as an OceanBioME `BiogeochemicalSediment` to pass to [`ECOSMO`](@ref)
as `sediment`. The pools start at zero; set them with `set!(sediment.fields.sed1, ...)`, in
mgC m⁻² like the Fortran.

With the default `bottom_stress = 0` the bottom stress never reaches `crBotStr`: detritus, opal and
calcite settle and nothing is resuspended.
"""
ecosmo_sediment(grid; parameters = EcosmoParameters(eltype(grid)), bottom_stress = 0) =
    BiogeochemicalSediment(grid, EcosmoSediment(parameters, convert_stress(eltype(grid), bottom_stress)))

convert_stress(FT, τ::Number) = convert(FT, τ)
convert_stress(FT, τ::QuadraticBottomStress) = QuadraticBottomStress{FT}(τ.drag_coefficient, τ.reference_density)

required_sediment_fields(::EcosmoSediment) = ECOSMO_SEDIMENT_FIELDS
required_tracers(s::EcosmoSediment) = (ECOSMO_BOTTOM_TRACKED_FIELDS..., stress_fields(s.bottom_stress)...)
sinking_fluxes(::EcosmoSediment) = ()
coupled_tracers(::EcosmoSediment) = ECOSMO_BOTTOM_COUPLED_TRACERS

stress_fields(::Number) = ()
stress_fields(::QuadraticBottomStress) = (:u, :v)

@inline bottom_stress(i, j, τ::Number, tracked) = τ

@inline function bottom_stress(i, j, τ::QuadraticBottomStress, tracked)
    u = @inbounds tracked.u[i, j, 1]
    v = @inbounds tracked.v[i, j, 1]
    return τ.reference_density * τ.drag_coefficient * (u * u + v * v)
end

@inline function column_bottom_fluxes(i, j, s::EcosmoSediment, fields, tracked)
    c = @inbounds (oxy = tracked.oxy[i, j, 1], no3 = tracked.no3[i, j, 1], det = tracked.det[i, j, 1],
                   opa = tracked.opa[i, j, 1], cal = tracked.cal[i, j, 1])
    pools = @inbounds NamedTuple{ECOSMO_SEDIMENT_FIELDS}(map(name -> fields[name][i, j, 1], ECOSMO_SEDIMENT_FIELDS))
    temp = @inbounds tracked.T[i, j, 1]
    om_cal = @inbounds tracked.Om_cal[i, j, 1]

    return bottom_fluxes(s.parameters, c, pools, temp, om_cal, bottom_stress(i, j, s.bottom_stress, tracked))
end

for name in ECOSMO_SEDIMENT_FIELDS
    @eval @inline (s::EcosmoSediment)(i, j, grid, ::Val{$(QuoteNode(name))}, clock, fields, tracked) =
        column_bottom_fluxes(i, j, s, fields, tracked).sediment.$name
end

for name in ECOSMO_BOTTOM_COUPLED_TRACERS
    @eval @inline (s::EcosmoSediment)(i, j, grid, ::Val{$(QuoteNode(name))}, clock, fields, tracked) =
        column_bottom_fluxes(i, j, s, fields, tracked).fluxes.$name
end

# OceanBioME adds a sediment flux to the bottom cell divided by `Δzᶜᶜᶠ`, the spacing between cell
# centres. A flux per unit area changes a concentration by flux / cell thickness — FABM's
# `cell_thickness` — and only that conserves mass on a stretched grid or a partial bottom cell.
function update_tendencies!(bgc, sediment::BiogeochemicalSediment{<:Any, <:EcosmoSediment}, model)
    grid = model.grid
    arch = architecture(grid)
    Gⁿ = model.timestepper.Gⁿ

    for name in coupled_tracers(sediment)
        args = (Val(name), sediment.biogeochemistry, sediment.fields, sediment.tracked_fields, model.clock)
        launch!(arch, grid, :xy, _add_bottom_flux!, Gⁿ[name], grid, sediment.bottom_indices, args)
    end

    return nothing
end

@kernel function _add_bottom_flux!(G, grid, bottom_indices, args)
    i, j = @index(Global, NTuple)

    @inbounds begin
        k = bottom_indices[i, j, 1]
        G[i, j, k] += sediment_tendencies(i, j, grid, args...) / Δzᶜᶜᶜ(i, j, k, grid)
    end
end

Base.summary(::EcosmoSediment{FT}) where {FT} = "ECOSMO sediment{$FT}"
