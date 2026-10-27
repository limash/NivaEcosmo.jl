"""
The NIVA adaptation of the ECOSMO marine biogeochemical model (`niva-ecosmo`, FABM/Fortran), ported
to Julia as an OceanBioME biogeochemistry for Oceananigans models.

The processes are pointwise functions shared by every driver — `interior_rates`, `bottom_fluxes`,
`oxygen_surface_flux`, `carbonate_diagnostics`, `carbon_dioxide_surface_flux`, `day_length`,
`gotm_light` — computing in the float type of the parameters and constants they are given.
[`ECOSMO`](@ref) assembles them into a model for Oceananigans, in the grid's float type.
"""
module NivaEcosmo

export ECOSMO,
    EcosmoParameters,
    GOTMLight,
    ecosmo_sediment,
    QuadraticBottomStress,
    ecosmo_surface_fluxes,
    fabm_configuration,
    ClipToMinimum,
    ECOSMO_TRACERS

using Adapt: Adapt, adapt
using KernelAbstractions: @kernel, @index
using YAML

using Oceananigans: Center, Face, Field, CenterField, ZFaceField, set!, fields
using Oceananigans.Architectures: architecture
using Oceananigans.BoundaryConditions: FluxBoundaryCondition, getbc
using Oceananigans.Fields: AbstractField, ZeroField
using Oceananigans.Grids: AbstractGrid, LatitudeLongitudeGrid, znode, φnode, xnode, ynode
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, immersed_cell
using Oceananigans.Operators: Δzᶜᶜᶜ
using Oceananigans.Utils: launch!
using Oceananigans.Biogeochemistry: AbstractBiogeochemistry

import Oceananigans.Biogeochemistry: required_biogeochemical_tracers,
                                     biogeochemical_auxiliary_fields,
                                     biogeochemical_drift_velocity,
                                     update_biogeochemical_state!,
                                     update_tendencies!

using OceanBioME: Biogeochemistry
using OceanBioME.Light: PAR_field
using OceanBioME.Models: teos10_polynomial_approximation
using OceanBioME.Sediments: BiogeochemicalSediment, AbstractSedimentBiogeochemistry, sediment_tendencies

import OceanBioME.Sediments: required_sediment_fields, required_tracers, sinking_fluxes, coupled_tracers

include("constants.jl")
include("parameters.jl")
include("ecosmo.jl")
include("carbonate.jl")
include("daylength.jl")
include("light.jl")
include("biogeochemistry.jl")
include("sediment.jl")
include("boundary_conditions.jl")
include("fabm_yaml.jl")

end # module
