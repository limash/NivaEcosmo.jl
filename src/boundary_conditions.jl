# The air-sea gas exchange of ECOSMO's `do_surface` (oxygen) and the carbonate module's `do_surface`
# (CO₂) as Oceananigans flux boundary conditions. FABM's surface fluxes are positive into the
# water; Oceananigans' are positive upward, hence the sign change.

"""
    ecosmo_surface_fluxes(FT = Float64; wind_speed, atmospheric_pCO2 = 420, pH_scale = 1,
                          density = teos10_polynomial_approximation)

Top boundary conditions for `oxy` and `dic`: `(oxy = FluxBoundaryCondition(...), dic = ...)`,
computing in `FT`, which must be the grid's float type.

`wind_speed` (m s⁻¹, at 10 m) is a number, a 2D field or a function `f(x, y, t)`.
`atmospheric_pCO2` is the mole fraction of CO₂ in air in ppm. `pH_scale` and `density` must be
the ones given to [`ECOSMO`](@ref). The CO₂ flux reads the `pCO2` auxiliary field, so it needs an
ECOSMO model.

Independent of the model object, so that a host can assemble its boundary conditions before it
builds the model, as Oceananigans requires.
"""
function ecosmo_surface_fluxes(FT::DataType = Float64; wind_speed, atmospheric_pCO2 = 420, pH_scale = 1,
                               density = teos10_polynomial_approximation)
    constants = ecosmo_constants(FT)
    wind_speed = convert_wind(FT, wind_speed)
    oxygen = FluxBoundaryCondition(oxygen_exchange; discrete_form = true, parameters = (; constants, wind_speed))

    carbon_parameters = (; constants, wind_speed, atmospheric_pCO2 = convert(FT, atmospheric_pCO2),
                         pH_scale = Int(pH_scale), density)
    carbon_dioxide = FluxBoundaryCondition(carbon_dioxide_exchange; discrete_form = true,
                                           parameters = carbon_parameters)

    return (oxy = oxygen, dic = carbon_dioxide)
end

convert_wind(FT, wind::Number) = convert(FT, wind)
convert_wind(FT, wind) = wind

@inline surface_value(value::Number, i, j, grid, clock) = value
@inline surface_value(field::AbstractField, i, j, grid, clock) = @inbounds field[i, j, 1]

@inline function surface_value(f::Function, i, j, grid, clock)
    x = xnode(i, j, grid.Nz, grid, Center(), Center(), Center())
    y = ynode(i, j, grid.Nz, grid, Center(), Center(), Center())
    return convert(eltype(grid), f(x, y, clock.time))  # whatever type `f` returns
end

@inline function oxygen_exchange(i, j, grid, clock, fields, p)
    k = grid.Nz
    wind = surface_value(p.wind_speed, i, j, grid, clock)
    flux = @inbounds oxygen_surface_flux(p.constants, fields.oxy[i, j, k], fields.T[i, j, k], fields.S[i, j, k], wind)
    return -flux
end

@inline function carbon_dioxide_exchange(i, j, grid, clock, fields, p)
    k = grid.Nz
    wind = surface_value(p.wind_speed, i, j, grid, clock)
    pressure = -znode(i, j, k, grid, Center(), Center(), Center())

    @inbounds begin
        T = fields.T[i, j, k]
        S = fields.S[i, j, k]
        flux = carbon_dioxide_surface_flux(p.constants, fields.dic[i, j, k], fields.alk[i, j, k], T, S,
                                           p.density(T, S, pressure / 10), pressure, wind,
                                           p.atmospheric_pCO2, fields.pCO2[i, j, k], p.pH_scale)
    end

    return -flux
end
