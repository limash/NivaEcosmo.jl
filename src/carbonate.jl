# The carbonate system of `niva_ecosmo_carbonate.F90` (the ERSEM carbonate module as adapted for
# ECOSMO), as pointwise functions of one cell. Only the branch the model uses is ported: `CO2DYN`
# always calls `CO2CLC` with `ICALC = 1` (DIC and alkalinity given) and boron on. The fits'
# coefficients are in `ECOSMO_CONSTANTS`.

# a / T + b + c ln T, the form of most of the fits below.
@inline temperature_fit((a, b, c), inverse_T, log_T) = a * inverse_T + b + c * log_T

# The factor by which pressure `P` (bar) changes an equilibrium constant at `T` (°C), from the
# changes in partial molar volume and compressibility (Millero 1995); `RT` in cm³ bar mol⁻¹.
@inline function pressure_correction(fit, T, P, RT)
    ΔV = evalpoly(T, fit.volume)
    Δκ = evalpoly(T, fit.compressibility)
    return exp((-ΔV + Δκ * P / 2) * P / RT)
end

# pK of carbonic acid, `fit` being (a0, a1, a2, b1, b2, b3, c1, c2, d) as in `ECOSMO_CONSTANTS`.
@inline function carbonic_acid_pK((a0, a1, a2, b1, b2, b3, c1, c2, d), inverse_T, log_T, sqrtS, S)
    return a0 + a1 * inverse_T + a2 * log_T + b1 * sqrtS + b2 * S + b3 * S^2 + (c1 * sqrtS + c2 * S) * inverse_T +
           d * sqrtS * log_T
end

# log10 of a calcium carbonate solubility product (Mucci 1983), `TK` in K.
@inline function mucci_log10_solubility(fit, TK, sqrtS, S)
    a0, a1, a2, a3 = fit.temperature
    b0, b1, b2 = fit.root_salinity
    c0, c1 = fit.salinity
    return a0 + a1 * TK + a2 / TK + a3 * log10(TK) + (b0 + b1 * TK + b2 / TK) * sqrtS + c0 * S + c1 * S * sqrtS
end

"""
    co2_constants(constants, P, T, S, pH_scale)

`CO2SET`: the equilibrium constants `(k0co2, k1co2, k2co2, kb)` at pressure `P` (bar), temperature
`T` (°C) and salinity `S`. `pH_scale` is the `pHscale` parameter: `1` total scale (Millero 2010),
`0` seawater scale (Millero 2010), `-1` the old seawater-scale formulation (Millero 1995).
"""
@inline function co2_constants(constants, P, T, S, pH_scale)
    TK = T + constants.zero_celsius
    inverse_T = 1 / TK
    log_T = log(TK)
    RT = constants.gas_constant * TK
    sqrtS = sqrt(S)

    # ionic strength, total sulfate and total fluoride
    a, b = constants.ionic_strength
    I = a * S / (1000 - b * S)
    chlorinity = S * constants.chlorinity_per_salinity
    sulfate = constants.sulfate_per_chlorinity * chlorinity
    fluoride = constants.fluoride_per_chlorinity * chlorinity

    # ks = [H][SO4]/[HSO4], free scale (Dickson 1990)
    bisulfate = constants.bisulfate
    ks = exp(temperature_fit(bisulfate.pure, inverse_T, log_T) +
             temperature_fit(bisulfate.root_ionic, inverse_T, log_T) * sqrt(I) +
             temperature_fit(bisulfate.ionic, inverse_T, log_T) * I +
             (bisulfate.ionic_1p5 * I * sqrt(I) + bisulfate.ionic_2 * I^2) * inverse_T +
             log(1 - bisulfate.salinity * S))

    # kf = [H][F]/[HF], total scale (Perez and Fraga 1987)
    f1, f2, f3 = constants.fluoride
    kf = exp(f1 * inverse_T + f2 + f3 * sqrtS)

    total2free_surface = 1 / (1 + sulfate / ks)

    # pressure corrections of ks and kf, the latter on the free scale (from Mocsy)
    ks = ks * pressure_correction(constants.bisulfate_pressure, T, P, RT)
    total2free_depth = 1 / (1 + sulfate / ks)
    kf = kf * total2free_surface * pressure_correction(constants.fluoride_pressure, T, P, RT) / total2free_depth
    total2sws_depth = total2free_depth * (1 + sulfate / ks + fluoride / (kf * total2free_depth))

    # k0co2 = CO2/fCO2 (Weiss 1974), with its pressure correction
    T100 = TK / 100
    solubility = constants.co2_solubility
    k0co2 = exp(temperature_fit(solubility.temperature, 1 / T100, log(T100)) +
                S * evalpoly(T100, solubility.salinity)) *
            exp((constants.standard_pressure - P) * constants.co2_partial_molar_volume / RT)

    # kb = [H][BO2]/[HBO2] (Millero 1995, from Dickson 1990)
    boric = constants.boric_acid
    kb = exp(evalpoly(sqrtS, boric.inverse_T) * inverse_T + evalpoly(sqrtS, boric.constant) +
             evalpoly(sqrtS, boric.log_T) * log_T + boric.T * sqrtS * TK)

    # k1co2 = [H][HCO3]/[H2CO3], k2co2 = [H][CO3]/[HCO3], on the pH scale chosen; kb is on the
    # total scale and moved to the seawater scale when that is the one
    scale = pH_scale + 2
    k1co2 = exp10(-carbonic_acid_pK(constants.carbonic_acid_1[scale], inverse_T, log_T, sqrtS, S))
    k2co2 = exp10(-carbonic_acid_pK(constants.carbonic_acid_2[scale], inverse_T, log_T, sqrtS, S))
    kb = kb * ifelse(pH_scale == 0, total2sws_depth, one(kb))

    # pressure corrections of k1co2, k2co2 and kb (Millero 1995)
    k1co2 = k1co2 * pressure_correction(constants.carbonic_acid_1_pressure, T, P, RT)
    k2co2 = k2co2 * pressure_correction(constants.carbonic_acid_2_pressure, T, P, RT)
    kb = kb * pressure_correction(constants.boric_acid_pressure, T, P, RT)

    return (; k0co2, k1co2, k2co2, kb)
end

# HALTAFALL gives up after this many passes; `CO2CLC`'s `C_CHECK`.
const CO2_MAXIMUM_ITERATIONS = 25

"""
    carbonate_system(constants, T, S, P, Ctot, TA, pH_scale)

`CO2DYN`: the carbonate speciation from DIC `Ctot` and alkalinity `TA` (both mol kg⁻¹) at
temperature `T` (°C), salinity `S` and pressure `P` (bar), by the HALTAFALL iteration of `CO2CLC`.

Returns `(; pH, pCO2, H2CO3, HCO3, CO3, k0co2, success)`, `pCO2` in atm and the species in
mol kg⁻¹. `success` is `false` when the iteration did not converge within
`CO2_MAXIMUM_ITERATIONS`; the other values are then meaningless, and the Fortran falls back on
the previous diagnostics.
"""
@inline function carbonate_system(constants, T, S, P, Ctot, TA, pH_scale)
    boron = constants.boron_per_salinity * S
    k = co2_constants(constants, P, max(T, 0), S, pH_scale)
    (; k1co2, k2co2, kb) = k

    # Each pass corrects carbonate and borate to their totals, then steps ln[H⁺] towards the
    # alkalinity: by a step that halves once the alkalinity has been bracketed, then by the secant
    # through the bracket.
    tolerance = constants.iteration_tolerance
    carbonate = TA / 10
    borate = boron
    hydrogen = constants.initial_hydrogen_ion

    phase = 1  # 1 at the start, ±2 stepping up or down, 3 once bracketed
    step_size = 2 * one(hydrogen)
    factor = one(hydrogen)
    x_above = y_above = x_below = y_below = zero(hydrogen)

    converged = false
    iterations = 0

    while !converged
        iterations += 1
        iterations >= CO2_MAXIMUM_ITERATIONS && break
        converged = true

        Y = carbonate * (1 + hydrogen / k2co2 + hydrogen^2 / (k1co2 * k2co2))
        if abs(Y - Ctot) > tolerance * Ctot
            carbonate = carbonate * Ctot / Y
            converged = false
        end

        Y = borate * (1 + hydrogen / kb)
        if abs(Y - boron) > tolerance * boron
            borate = borate * boron / Y
            converged = false
        end

        # Alkalinity is equivalent to -(total H⁺), so the sign of `direction` is opposite to the
        # usual one.
        Y = carbonate * (2 + hydrogen / k2co2) + borate
        if abs(Y - TA) > tolerance * TA
            converged = false
            x = log(hydrogen)
            direction = sign(Y - TA)
            if direction >= 0
                x_above, y_above = x, Y
            else
                x_below, y_below = x, Y
            end

            previous_phase = phase
            if previous_phase == 1
                phase = ifelse(direction > 0, 2, -2)
            elseif abs(previous_phase) == 2 && previous_phase * direction < 0
                factor = one(factor) / 2
                phase = 3
            end

            if phase == 3 && step_size < constants.secant_threshold
                x = x_above + (x_below - x_above) / (y_below - y_above) * (TA - y_above)
            else
                step_size = step_size * factor
                x = x + step_size * direction
            end
            hydrogen = exp(x)
        end
    end

    HCO3 = carbonate * hydrogen / k2co2
    H2CO3 = HCO3 * hydrogen / k1co2
    pCO2 = H2CO3 / k.k0co2
    pH = -log10(hydrogen)

    return (; pH, pCO2, H2CO3, HCO3, CO3 = carbonate, k0co2 = k.k0co2, success = converged)
end

"""
    calcite_aragonite_saturation(constants, T, S, P, CO3)

`CaCO3_Saturation`: the saturation states `(Om_cal, Om_arg)` of calcite and aragonite at
temperature `T` (°C), salinity `S`, pressure `P` (bar) and carbonate ion concentration `CO3`
(mol kg⁻¹), after Zeebe and Wolf-Gladrow (2001) and Mucci (1983) with Millero's (1995) pressure
correction.
"""
@inline function calcite_aragonite_saturation(constants, T, S, P, CO3)
    TK = T + constants.zero_celsius
    RT = constants.gas_constant * TK
    sqrtS = sqrt(S)

    Kspc = exp10(mucci_log10_solubility(constants.calcite_solubility, TK, sqrtS, S)) *
           pressure_correction(constants.calcite_pressure, T, P, RT)
    Kspa = exp10(mucci_log10_solubility(constants.aragonite_solubility, TK, sqrtS, S)) *
           pressure_correction(constants.aragonite_pressure, T, P, RT)

    ion_product = CO3 * constants.calcium
    return (ion_product / Kspc, ion_product / Kspa)
end

"""
    carbonate_diagnostics(constants, dic, alk, temp, salt, density, pressure, previous, pH_scale)

The carbonate module's `do` at one cell: `(; pH, pCO2, CarbA, BiCarb, Carb, Om_cal, Om_arg)`, with
`pCO2` in µatm and the species in mmol m⁻³, as FABM reports them.

`dic` and `alk` are in mmol m⁻³, `density` in kg m⁻³ and `pressure` in dbar. When the iteration
does not converge, the Fortran keeps its previous diagnostics and recomputes only the saturation
states, from the previous carbonate ion concentration; `previous` supplies those
(`pH`, `pCO2`, `CarbA`, `BiCarb`, `Carb`), and the result is the same in that case.
"""
@inline function carbonate_diagnostics(constants, dic, alk, temp, salt, density, pressure, previous, pH_scale)
    S = max(salt, 0)
    P = pressure / 10  # dbar → bar
    to_mol_per_kg = 1 / (1000 * density)  # from mmol m⁻³

    cs = carbonate_system(constants, temp, S, P, dic * to_mol_per_kg, alk * to_mol_per_kg, pH_scale)

    if cs.success
        pH = cs.pH
        pCO2 = cs.pCO2 * 1_000_000  # atm → µatm
        CarbA = cs.H2CO3 / to_mol_per_kg
        BiCarb = cs.HCO3 / to_mol_per_kg
        Carb = cs.CO3 / to_mol_per_kg
    else
        (; pH, pCO2, CarbA, BiCarb, Carb) = previous
    end

    Om_cal, Om_arg = calcite_aragonite_saturation(constants, temp, S, P, Carb * to_mol_per_kg)

    return (; pH, pCO2, CarbA, BiCarb, Carb, Om_cal, Om_arg)
end

"""
    carbon_dioxide_surface_flux(constants, dic, alk, temp, salt, density, pressure, wind, atmospheric_pCO2,
                                previous_pCO2, pH_scale)

The carbonate module's `do_surface`: the air-sea CO₂ flux in mmol C m⁻² s⁻¹, positive into the
water, with the Wanninkhof (2014) transfer velocity.

`atmospheric_pCO2` is the mole fraction of CO₂ in air in ppm. `previous_pCO2` (µatm) is used when
the carbonate iteration does not converge, as in the Fortran.
"""
@inline function carbon_dioxide_surface_flux(constants, dic, alk, temp, salt, density, pressure, wind,
                                             atmospheric_pCO2, previous_pCO2, pH_scale)
    to_mol_per_kg = 1 / (1000 * density)  # from mmol m⁻³
    cs = carbonate_system(constants, temp, max(salt, 0), pressure / 10, dic * to_mol_per_kg, alk * to_mol_per_kg,
                          pH_scale)
    pCO2 = ifelse(cs.success, cs.pCO2, previous_pCO2 / 1_000_000)  # atm

    velocity = transfer_velocity(constants, constants.schmidt_CO2, temp, wind)
    return velocity * cs.k0co2 * (atmospheric_pCO2 / 1_000_000 - pCO2) / to_mol_per_kg
end
