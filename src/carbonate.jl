# The carbonate system of `niva_ecosmo_carbonate.F90` (the ERSEM carbonate module as adapted for
# ECOSMO), as pointwise functions of one cell. Only the branch the model uses is ported: `CO2DYN`
# always calls `CO2CLC` with `ICALC = 1` (DIC and alkalinity given) and boron on.
#
# Two literals in `CO2SET` have no kind suffix and are therefore *single* precision in the Fortran
# — `771.54` and `5.7789` — and are written as `Float32` constants here so they round the same.

const RGAS_CO2 = 83.131  # cm³ bar mol⁻¹ K⁻¹

"""
    co2_constants(P, T, S, pH_scale)

`CO2SET`: the equilibrium constants `(k0co2, k1co2, k2co2, kb)` at pressure `P` (bar), temperature
`T` (°C) and salinity `S`. `pH_scale` is the `pHscale` parameter: `1` total scale (Millero 2010),
`0` seawater scale (Millero 2010), `-1` the old seawater-scale formulation (Millero 1995).
"""
@inline function co2_constants(P, T, S, pH_scale)
    TK = T + 273.15
    invtk = 1.0 / TK
    dlogTK = log(TK)
    TK100 = TK / 100.0
    S2 = S * S
    sqrtS = sqrt(S)
    S15 = S^1.5

    # ionic strength, total sulfate and total fluoride
    is = 19.924 * S / (1000.0 - 1.005 * S)
    sqrtis = sqrt(is)
    cl = S / 1.80655
    st = 0.14 * cl / 96.065

    # ks = [H][SO4]/[HSO4], free scale (Dickson 1990)
    ks = exp(-4276.1 * invtk + 141.328 - 23.093 * dlogTK +
             (-13856.0 * invtk + 324.57 - 47.986 * dlogTK) * sqrtis +
             (35474.0 * invtk - Float64(771.54f0) + 114.723 * dlogTK) * is -
             2698.0 * invtk * is^1.5 + 1776.0 * invtk * (is * is) +
             log(1.0 - 0.001005 * S))

    ft = 0.000067 * cl / 19.9984

    # kf = [H][F]/[HF], total scale (Perez and Fraga 1987)
    kf = exp(874.0 * invtk - 9.68 + 0.111 * sqrtS)

    total2free_surface = 1.0 / (1.0 + st / ks)
    free2sws_surface = 1.0 + st / ks + ft / kf
    total2sws_surface = total2free_surface * free2sws_surface

    # pressure correction of ks (from Mocsy)
    delta = -18.03 + 0.0466 * T + 0.000316 * (T * T)
    kappa = -4.53 + 0.00009 * T
    ks = ks * exp((-delta + 0.5 * kappa * P) * P / (RGAS_CO2 * TK))
    total2free_depth = 1.0 / (1.0 + st / ks)

    # pressure correction of kf, which needs kf on the free scale
    delta = -9.78 - 0.009 * T - 0.000942 * (T * T)
    kappa = -3.91 + 0.000054 * T
    kf = kf * total2free_surface * exp((-delta + 0.5 * kappa * P) * P / (RGAS_CO2 * TK)) / total2free_depth
    free2sws_depth = 1.0 + st / ks + ft / (kf * total2free_depth)
    total2sws_depth = total2free_depth * free2sws_depth

    # k0co2 = CO2/fCO2 (Weiss 1974), with its pressure correction
    k0co2 = exp(93.4517 / TK100 - 60.2409 + 23.3585 * log(TK100) +
                S * (0.023517 - 0.023656 * TK100 + 0.0047036 * (TK100 * TK100)))
    k0co2 = k0co2 * exp(((1.01325 - P) * 32.3) / (RGAS_CO2 * TK))

    # kb = [H][BO2]/[HBO2] (Millero 1995, from Dickson 1990)
    kb = exp((-8966.9 - 2890.53 * sqrtS - 77.942 * S + 1.728 * S15 - 0.0996 * S2) / TK +
             (148.0248 + 137.1942 * sqrtS + 1.62142 * S) +
             (-24.4344 - 25.085 * sqrtS - 0.2474 * S) * dlogTK +
             0.053105 * sqrtS * TK)

    # k1co2 = [H][HCO3]/[H2CO3], k2co2 = [H][CO3]/[HCO3]
    if pH_scale == -1
        k1co2 = 10.0^(-1.0 * (3670.7 / TK - 62.008 + 9.7944 * dlogTK - 0.0118 * S + 0.000116 * S2))
        k2co2 = 10.0^(-1.0 * (1394.7 / TK + 4.777 - 0.0184 * S + 0.000118 * S2))
    elseif pH_scale == 0
        pk1co2 = -126.34048 + 6320.813 * invtk + 19.568224 * dlogTK +
                 (13.4038 * sqrtS + 0.03206 * S - 0.00005242 * (S * S)) +
                 (-530.659 * sqrtS - 5.8210 * S) * invtk +
                 (-2.0664 * sqrtS) * dlogTK
        k1co2 = 10.0^(-pk1co2)
        pk2co2 = -90.18333 + 5143.692 * invtk + 14.613358 * dlogTK +
                 (21.3728 * sqrtS + 0.1218 * S - 0.0003688 * (S * S)) +
                 (-788.289 * sqrtS - 19.189 * S) * invtk +
                 (-3.374 * sqrtS) * dlogTK
        k2co2 = 10.0^(-pk2co2)
        kb = kb * total2sws_depth
    else
        pk1co2 = -126.34048 + 6320.813 * invtk + 19.568224 * dlogTK +
                 (13.4051 * sqrtS + 0.03185 * S - 0.00005218 * (S * S)) +
                 (-531.095 * sqrtS - Float64(5.7789f0) * S) * invtk +
                 (-2.0663 * sqrtS) * dlogTK
        k1co2 = 10.0^(-pk1co2)
        pk2co2 = -90.18333 + 5143.692 * invtk + 14.613358 * dlogTK +
                 (21.5724 * sqrtS + 0.1212 * S - 0.0003714 * (S * S)) +
                 (-798.292 * sqrtS - 18.951 * S) * invtk +
                 (-3.403 * sqrtS) * dlogTK
        k2co2 = 10.0^(-pk2co2)
    end

    # pressure corrections of k1co2, k2co2 and kb (Millero 1995)
    delta = -25.5 + 0.1271 * T
    kappa = (-3.08 + 0.0877 * T) / 1000.0
    k1co2 = k1co2 * exp((-delta + 0.5 * kappa * P) * P / (RGAS_CO2 * TK))

    delta = -15.82 - 0.0219 * T
    kappa = (1.13 - 0.1475 * T) / 1000.0
    k2co2 = k2co2 * exp((-delta + 0.5 * kappa * P) * P / (RGAS_CO2 * TK))

    delta = -29.48 + 0.1622 * T - 0.002608 * (T * T)
    kappa = -2.84 / 1000.0
    kb = kb * exp((-delta + 0.5 * kappa * P) * P / (RGAS_CO2 * TK))

    return (; k0co2, k1co2, k2co2, kb)
end

# HALTAFALL gives up after this many passes; `CO2CLC`'s `C_CHECK`.
const CO2_MAXIMUM_ITERATIONS = 25

"""
    carbonate_system(T, S, P, Ctot, TA, pH_scale)

`CO2DYN`: the carbonate speciation from DIC `Ctot` and alkalinity `TA` (both mol kg⁻¹) at
temperature `T` (°C), salinity `S` and pressure `P` (bar), by the HALTAFALL iteration of `CO2CLC`.

Returns `(; pH, pCO2, H2CO3, HCO3, CO3, k0co2, success)`, `pCO2` in atm and the species in
mol kg⁻¹. `success` is `false` when the iteration did not converge within
`CO2_MAXIMUM_ITERATIONS`; the other values are then meaningless, and the Fortran falls back on
the previous diagnostics.
"""
@inline function carbonate_system(T, S, P, Ctot, TA, pH_scale)
    Tmax = max(T, 0.0)
    Btot = 0.0004128 * S / 35.0

    k = co2_constants(P, Tmax, S, pH_scale)
    k1co2 = k.k1co2
    k2co2 = k.k2co2
    kb = k.kb

    # initial guesses and tolerances (PCO2 enters as 0, hence H2CO3 = 0)
    CO3 = TA / 10.0
    AHPLUS = 1.0e-8
    ALKB = Btot
    TOL1 = TA / 1.0e5
    TOL3 = Ctot / 1.0e5
    TOL4 = Btot / 1.0e5

    KARL = 1
    STEG = 2.0
    FAK = 1.0
    STEGBY = 0.4
    X1 = Y1 = X2 = Y2 = 0.0

    done = false
    counter = 0
    success = true

    while !done
        done = true
        counter += 1

        if counter >= CO2_MAXIMUM_ITERATIONS
            success = false
            break
        end

        # Ctot is fixed
        Y = CO3 * (1.0 + AHPLUS / k2co2 + AHPLUS * AHPLUS / (k1co2 * k2co2))
        if abs(Y - Ctot) > TOL3
            CO3 = CO3 * Ctot / Y
            done = false
        end

        Y = ALKB * (1.0 + AHPLUS / kb)
        if abs(Y - Btot) > TOL4
            ALKB = ALKB * Btot / Y
            done = false
        end

        # Alkalinity is equivalent to -(total H+), so the sign of W is opposite to the usual one.
        Y = CO3 * (2.0 + AHPLUS / k2co2) + ALKB
        if abs(Y - TA) > TOL1
            done = false
            X = log(AHPLUS)
            W = copysign(1.0, Y - TA)
            if W >= 0.0
                X1 = X
                Y1 = Y
            else
                X2 = X
                Y2 = Y
            end
            LQ = KARL
            if LQ == 1
                KARL = 2 * unsafe_trunc(Int, W)  # NINT(W), W being ±1
            elseif abs(LQ) == 2 && LQ * W < 0.0
                FAK = 0.5
                KARL = 3
            end
            if KARL == 3 && STEG < STEGBY
                W = (X2 - X1) / (Y2 - Y1)
                X = X1 + W * (TA - Y1)
            else
                STEG = STEG * FAK
                X = X + STEG * W
            end
            AHPLUS = exp(X)
        end
    end

    HCO3 = CO3 * AHPLUS / k2co2
    H2CO3 = HCO3 * AHPLUS / k1co2
    pCO2 = H2CO3 / k.k0co2
    pH = -log10(AHPLUS)

    return (; pH, pCO2, H2CO3, HCO3, CO3, k0co2 = k.k0co2, success)
end

"""
    calcite_aragonite_saturation(Tc, S, Pr, CO3)

`CaCO3_Saturation`: the saturation states `(Om_cal, Om_arg)` of calcite and aragonite at
temperature `Tc` (°C), salinity `S`, pressure `Pr` (Pa) and carbonate ion concentration `CO3`
(mol kg⁻¹), after Zeebe and Wolf-Gladrow (2001) and Mucci (1983) with Millero's (1995) pressure
correction.
"""
@inline function calcite_aragonite_saturation(Tc, S, Pr, CO3)
    Tk = Tc + 273.15
    Ca = 0.01028
    R = 83.131
    P = Pr * 1.0e-5

    # calcite
    tmp1 = -171.9065 - (0.077993 * Tk) + (2839.319 / Tk) + 71.595 * log10(Tk)
    tmp2 = (-0.77712 + (0.0028426 * Tk) + (178.34 / Tk)) * sqrt(S)
    tmp3 = -(0.07711 * S) + (0.0041249 * (S^1.5))
    Kspc = 10.0^(tmp1 + tmp2 + tmp3)

    dV = -48.76 + 0.5304 * Tc
    dK = -11.76 / 1.0e3 + (0.3692 / 1.0e3) * Tc
    tmp1 = -(dV / (R * Tk)) * P + (0.5 * dK / (R * Tk)) * P * P
    Kspc = Kspc * exp(tmp1)

    # aragonite
    tmp1 = -171.945 - 0.077993 * Tk + 2903.293 / Tk + 71.595 * log10(Tk)
    tmp2 = (-0.068393 + 0.0017276 * Tk + 88.135 / Tk) * sqrt(S)
    tmp3 = -0.10018 * S + 0.0059415 * S^1.5
    Kspa = 10.0^(tmp1 + tmp2 + tmp3)

    dV = -46.0 + 0.530 * Tc
    dK = -11.76 / 1.0e3 + (0.3692 / 1.0e3) * Tc
    tmp1 = -(dV / (R * Tk)) * P + (0.5 * dK / (R * Tk)) * P * P
    Kspa = Kspa * exp(tmp1)

    return ((CO3 * Ca) / Kspc, (CO3 * Ca) / Kspa)
end

"""
    carbonate_diagnostics(dic, alk, temp, salt, density, pressure, previous, pH_scale)

The carbonate module's `do` at one cell: `(; pH, pCO2, CarbA, BiCarb, Carb, Om_cal, Om_arg)`, with
`pCO2` in µatm and the species in mmol m⁻³, as FABM reports them.

`dic` and `alk` are in mmol m⁻³, `density` in kg m⁻³ and `pressure` in dbar. When the iteration
does not converge, the Fortran keeps its previous diagnostics and recomputes only the saturation
states, from the previous carbonate ion concentration; `previous` supplies those
(`pH`, `pCO2`, `CarbA`, `BiCarb`, `Carb`), and the result is the same in that case.
"""
@inline function carbonate_diagnostics(dic, alk, temp, salt, density, pressure, previous, pH_scale)
    X1X = max(salt, 0.0)
    TA = alk / 1.0e3 / density      # mmol m⁻³ → mol kg⁻¹
    Ctot = dic / 1.0e3 / density

    cs = carbonate_system(temp, X1X, pressure * 0.1, Ctot, TA, pH_scale)  # dbar → bar

    if cs.success
        pH = cs.pH
        pCO2 = cs.pCO2 * 1.0e6
        CarbA = cs.H2CO3 * 1.0e3 * density
        BiCarb = cs.HCO3 * 1.0e3 * density
        Carb = cs.CO3 * 1.0e3 * density
        CO3 = cs.CO3
    else
        pH = previous.pH
        pCO2 = previous.pCO2
        CarbA = previous.CarbA
        BiCarb = previous.BiCarb
        Carb = previous.Carb
        CO3 = previous.Carb / 1.0e3 / density
    end

    Om_cal, Om_arg = calcite_aragonite_saturation(temp, X1X, pressure * 1.0e4, CO3)  # dbar → Pa

    return (; pH, pCO2, CarbA, BiCarb, Carb, Om_cal, Om_arg)
end

"""
    carbon_dioxide_surface_flux(dic, alk, temp, salt, density, pressure, wind, atmospheric_pCO2,
                                previous_pCO2, pH_scale)

The carbonate module's `do_surface`: the air-sea CO₂ flux in mmol C m⁻² s⁻¹, positive into the
water, with the Wanninkhof (2014) transfer velocity.

`atmospheric_pCO2` is the mole fraction of CO₂ in air in ppm. `previous_pCO2` (µatm) is used when
the carbonate iteration does not converge, as in the Fortran.
"""
@inline function carbon_dioxide_surface_flux(dic, alk, temp, salt, density, pressure, wind,
                                             atmospheric_pCO2, previous_pCO2, pH_scale)
    S = max(salt, 0.0)
    wnd = max(wind, 0.0)

    TA = alk / 1.0e3 / density
    Ctot = dic / 1.0e3 / density

    cs = carbonate_system(temp, S, pressure * 0.1, Ctot, TA, pH_scale)
    PCO2 = cs.success ? cs.pCO2 : previous_pCO2 * 1.0e-6

    # Schmidt number for CO2 following Wanninkhof (2014)
    T = max(min(temp, 40.0), -2.0)
    sc = 2116.8 - 136.25 * T + 4.7353 * (T * T) - 0.092307 * (T * T * T) + 0.0007555 * T^4.0

    fwind = 0.251 * (wnd * wnd) * (sc / 660.0)^(-0.5)  # cm/hr
    fwind = fwind / 360000.0                            # m/s
    UPTAKE = fwind * cs.k0co2 * (atmospheric_pCO2 / 1.0e6 - PCO2)

    return UPTAKE * 1.0e3 * density
end
