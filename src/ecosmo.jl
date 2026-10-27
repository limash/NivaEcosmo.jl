# The ECOSMO processes as pointwise functions of one cell, ported from `niva_ecosmo.F90`.
#
# These are what both drivers call: the OceanBioME tendency kernels and the box model. Every number
# in them is a parameter, an entry of the constants (`ecosmo_constants`), or an exact integer, so
# they compute in the float type of the parameters and constants they are given.

"""
    interior_rates(p::EcosmoParameters, constants, c, temp, D, E0, om_cal)

The interior sources of every ECOSMO tracer and the diagnostics of `niva_ecosmo`'s `do`, at one point.

`constants` is [`ecosmo_constants`](@ref) in the parameters' float type. `c` holds the
concentrations by FABM name (`no3`, `nh4`, `pho`, `sil`, `oxy`, `fla`, `dia`, `flachl`, `diachl`,
`microzoo`, `mesozoo`, `det`, `opa`, `cal`, `dom`), `temp` is in °C, `D` the day length as a
fraction of the day, `E0` the scalar PAR in W m⁻² and `om_cal` the calcite saturation state.

Returns `(; rates, diagnostics)`, `rates` per second and including `dic` and `alk`, which the
Fortran adds to the carbonate module's state variables when `couple_co2` is on.
"""
@inline function interior_rates(p::EcosmoParameters, constants, c, temp, D, E0, om_cal)
    (; no3, nh4, pho, sil, oxy, fla, dia, flachl, diachl, microzoo, mesozoo, det, opa, cal, dom) = c
    (; redfield_N_C, redfield_P_C, redfield_Si_C, mol_per_gram_C, nitrate_oxygen_ratio) = constants
    ε = constants.uptake_regularization

    E0c = max(E0, 0)

    # Loss terms are switched off below a threshold biomass, to keep plankton from extinction.
    fla_loss = at_least(fla, constants.phytoplankton_loss_threshold)
    dia_loss = at_least(dia, constants.phytoplankton_loss_threshold)
    mic_loss = at_least(microzoo, constants.zooplankton_loss_threshold)
    mes_loss = at_least(mesozoo, constants.zooplankton_loss_threshold)

    # nutrient limitation factors
    up_ndrPs = p.rNO3Ps * p.rNH4Ps + p.rNH4Ps * no3 + p.rNO3Ps * nh4
    up_ndrPl = p.rNO3Pl * p.rNH4Pl + p.rNH4Pl * no3 + p.rNO3Pl * nh4
    up_nh4Ps = p.rNO3Ps * nh4 / up_ndrPs
    up_nh4Pl = p.rNO3Pl * nh4 / up_ndrPl
    up_no3Ps = p.rNH4Ps * no3 / up_ndrPs
    up_no3Pl = p.rNH4Pl * no3 / up_ndrPl
    up_nPs = up_nh4Ps + up_no3Ps
    up_nPl = up_nh4Pl + up_no3Pl
    up_phoPs = pho / (p.rPO4Ps + pho)
    up_phoPl = pho / (p.rPO4Pl + pho)
    t_silPl = max(sil - p.minSiPl, 0)
    up_silPl = t_silPl / (p.rSiPl + t_silPl)
    NutlimPs = min(up_nPs, up_phoPs)
    NutlimPl = min(up_nPl, up_phoPl, up_silPl)

    # maximum phytoplankton growth rates, corrected for temperature and fractional day length
    D₀ = constants.reference_day_length
    muPsc = q10_factor(p.q10Ps, temp) * (p.gammaDPs + D₀) / (p.gammaDPs + D) * p.muPs
    muPlc = q10_factor(p.q10Pl, temp) * (p.gammaDPl + D₀) / (p.gammaDPl + D) * p.muPl

    # further correction for light limitation
    ChlCPs = clamp(flachl / max(fla, constants.minimum_carbon), p.MINchl2cPs, p.MAXchl2cPs)
    ChlCPl = clamp(diachl / max(dia, constants.minimum_carbon), p.MINchl2cPl, p.MAXchl2cPl)

    if E0 > constants.minimum_light
        muPsc = muPsc * (1 - exp(-p.alfaPs * E0 * ChlCPs / muPsc)) * exp(-p.betaPs * E0 * ChlCPs / muPsc)
        muPlc = muPlc * (1 - exp(-p.alfaPl * E0 * ChlCPl / muPlc)) * exp(-p.betaPl * E0 * ChlCPl / muPlc)

        # ratio of chlorophyll-a synthesis to carbon fixation rates
        chl2c_fla = clamp(p.MAXchl2cPs * (muPsc * fla) / (p.alfaPs * E0 * flachl), p.MINchl2cPs, p.MAXchl2cPs)
        chl2c_dia = clamp(p.MAXchl2cPl * (muPlc * dia) / (p.alfaPl * E0 * diachl), p.MINchl2cPl, p.MAXchl2cPl)
    else
        # The E0 → 0 limit of the expressions above (L'Hôpital).
        muPsc = zero(muPsc)
        muPlc = zero(muPlc)
        chl2c_fla = p.MAXchl2cPs
        chl2c_dia = p.MAXchl2cPl
    end

    # further correction for nutrient limitation
    muPsc = muPsc * NutlimPs
    muPlc = muPlc * NutlimPl

    # primary productivity
    ProdPs = muPsc * fla
    ProdPl = muPlc * dia
    Prod = ProdPs + ProdPl

    # phytoplankton mortality, enhanced under nutrient limitation
    NutfacmPs = p.amortPs / (1 + NutlimPs * (p.amortPs - 1))
    NutfacmPl = p.amortPl / (1 + NutlimPl * (p.amortPl - 1))
    mort1Ps = NutfacmPs * p.m1Ps * fla_loss * fla
    mort1Pl = NutfacmPl * p.m1Pl * dia_loss * dia
    mort2Ps = NutfacmPs * p.m2Ps * fla * fla_loss * fla
    mort2Pl = NutfacmPl * p.m2Pl * dia * dia_loss * dia
    mortPs = mort1Ps + mort2Ps
    mortPl = mort1Pl + mort2Pl

    # grazing/feeding fluxes
    Fs = p.prefZsPs * fla + p.prefZsPl * dia + p.prefZsZs * microzoo + p.prefZsD * det
    Fl = p.prefZlPs * fla + p.prefZlPl * dia + p.prefZlZs * microzoo + p.prefZlZl * mesozoo + p.prefZlD * det

    etZs = q10_factor(p.q10Zs, temp)
    etZl = q10_factor(p.q10Zl, temp)

    ZsonPs = fla_loss * p.GrZsPs * etZs * p.prefZsPs * fla / (p.RgZs + Fs) * microzoo
    ZsonPl = dia_loss * p.GrZsPl * etZs * p.prefZsPl * dia / (p.RgZs + Fs) * microzoo
    ZsonZs = mic_loss * p.GrZsZs * etZs * p.prefZsZs * microzoo / (p.RgZs + Fs) * microzoo
    ZsonD = p.GrZsD * etZs * p.prefZsD * det / (p.RgZs + Fs) * microzoo

    ZlonPs = fla_loss * p.GrZlPs * etZl * p.prefZlPs * fla / (p.RgZl + Fl) * mesozoo
    ZlonPl = dia_loss * p.GrZlPl * etZl * p.prefZlPl * dia / (p.RgZl + Fl) * mesozoo
    ZlonZs = mic_loss * p.GrZlZs * etZl * p.prefZlZs * microzoo / (p.RgZl + Fl) * mesozoo
    ZlonZl = mes_loss * p.GrZlZl * etZl * p.prefZlZl * mesozoo / (p.RgZl + Fl) * mesozoo
    ZlonD = p.GrZlD * etZl * p.prefZlD * det / (p.RgZl + Fl) * mesozoo

    f_O2pos, f_O2resp, f_denit, f_SO4resp = redox_switches(oxy, no3)

    # phytoplankton
    rhs_fla = ProdPs - mortPs - ZsonPs - ZlonPs
    rhs_dia = ProdPl - mortPl - ZsonPl - ZlonPl

    # phytoplankton chlorophyll-a
    rhs_flachl = ProdPs * chl2c_fla - (mortPs + ZsonPs + ZlonPs) * flachl / fla
    rhs_diachl = ProdPl * chl2c_dia - (mortPl + ZsonPl + ZlonPl) * diachl / dia

    # microzooplankton
    egesZs = (1 - p.gammaZsP) * (ZsonPs + ZsonPl + ZsonZs) + (1 - p.gammaZsD) * ZsonD
    excrZs = p.excZs * etZs * mic_loss * microzoo + p.fAexcZs * (p.gammaZsP * (ZsonPs + ZsonPl + ZsonZs) + p.gammaZsD * ZsonD)
    etmZs = q10_factor(p.q10mZs, temp)
    mortZs = p.mZs * etmZs * mic_loss * microzoo
    prodZs = ZsonPs + ZsonPl + ZsonZs + ZsonD - egesZs
    rhs_microzoo = prodZs - excrZs - mortZs - ZsonZs - ZlonZs

    # mesozooplankton, ingestion corrected for sloppy feeding
    inZlPs = (1 - p.fGslpZlPs) * ZlonPs
    inZlPl = (1 - p.fGslpZlPl) * ZlonPl
    inZlZs = (1 - p.fGslpZlZs) * ZlonZs
    inZlZl = (1 - p.fGslpZlZl) * ZlonZl
    inZlPZ = inZlPs + inZlPl + inZlZs + inZlZl
    inZlD = (1 - p.fGslpZlD) * ZlonD
    egesZl = (1 - p.gammaZlP) * inZlPZ + (1 - p.gammaZlD) * inZlD
    excrZl = p.excZl * etZl * mes_loss * mesozoo + p.fAexcZl * (p.gammaZlP * inZlPZ + p.gammaZlD * inZlD)
    etmZl = q10_factor(p.q10mZl, temp)
    mZlt = p.mZl + p.mvZl * E0c / (p.KE0mvZl + E0c)  # total mortality rate inc. visual predation
    mortZl = mZlt * etmZl * mes_loss * mesozoo
    prodZl = inZlPZ + inZlD - egesZl
    rhs_mesozoo = prodZl - excrZl - mortZl - ZlonZl

    # detritus
    dxxt = mortPs + mortPl + mortZs + mortZl + egesZs + egesZl
    dxxdom = p.frmort1Ps * mort1Ps + p.frmort1Pl * mort1Pl +
             p.frmort2Ps * mort2Ps + p.frmort2Pl * mort2Pl +
             p.frmortZs * mortZs + p.frmortZl * mortZl +
             p.fregesZs * egesZs + p.fregesZl * egesZl
    frem = p.reminD * (1 + 20 * (temp^2 / (13^2 + temp^2)))
    remdet = frem * det
    rhs_det = (dxxt - dxxdom) - ZsonD - ZlonD - remdet

    # labile dissolved organic matter
    slpdom = (p.fGslpZlPs * ZlonPs + p.fGslpZlPl * ZlonPl + p.fGslpZlZs * ZlonZs + p.fGslpZlZl * ZlonZl +
              p.fGslpZlD * ZlonD)
    excrt = excrZs + excrZl
    excrdom = p.fexcdomZs * excrZs + p.fexcdomZl * excrZl
    excrmin = excrt - excrdom
    fremDOM = 10 * frem
    remdom = fremDOM * dom
    rhs_dom = dxxdom + slpdom + excrdom - remdom

    # nitrate
    nitrif = nh4 * p.nitrifmax * exp(p.aTnitrif * temp) *  # T-dependence
             f_O2pos * oxy / (p.KO2nitrif + oxy) *          # O2-dependence
             1 / (1 + p.aE0nitrif * E0c^p.bE0nitrif)        # PAR-dependence
    remtot = remdet + remdom
    denit = f_denit * p.eta_denit * remtot
    rhs_nit = -(up_no3Ps + ε / 2) / (up_nPs + ε) * ProdPs -
              (up_no3Pl + ε / 2) / (up_nPl + ε) * ProdPl +
              nitrif -
              denit

    # ammonium
    rhs_amm = -(up_nh4Ps + ε / 2) / (up_nPs + ε) * ProdPs -
              (up_nh4Pl + ε / 2) / (up_nPl + ε) * ProdPl +
              excrmin +
              remtot -
              nitrif

    # phosphate
    rhs_pho = -Prod + excrmin + remtot

    # silicate: opal dissolution after Ridgwell et al. (2002)
    a, b = constants.opal_solubility
    sil_eq = exp10(a - b / (temp + constants.zero_celsius))
    u_opal = clamp((sil_eq - sil * redfield_Si_C * mol_per_gram_C) / sil_eq, 0, 1)
    dissolution = constants.opal_dissolution
    ndissO = dissolution.linear * (1 + temp / 15) * u_opal +
             dissolution.power * (1 + temp / 400)^37 * u_opal^dissolution.exponent
    dissO = p.dissO * ndissO * opa
    rhs_sil = -ProdPl + dissO

    # opal
    rhs_opa = mortPl + ZsonPl + ZlonPl - dissO

    # calcite: rain ratio within calcifiers, scaled to the small phytoplankton as a whole
    RainR = p.Rain0 * max(0, (om_cal - 1) / (om_cal - 1 + p.Kcalom))
    t = max(0, temp)
    calc_Ps = min(1 - up_phoPs, up_nPs) * (t / (2 + t))  # calcifying fraction of Ps
    RainR = max(RainR * calc_Ps, constants.minimum_rain_ratio)
    calc = RainR * (((1 - p.frmort1Ps) * mort1Ps + (1 - p.frmort2Ps) * mort2Ps) +
                    ZsonPs * (1 - p.fdissCZs) +
                    ZlonPs * (1 - p.fdissCZl))
    diss = p.dissCmax * max(1 - om_cal, 0)^p.ndissC * cal
    rhs_cal = calc - diss

    # oxygen
    rhs_oxy = ((up_nh4Ps + nitrate_oxygen_ratio * up_no3Ps + ε) / (up_nPs + ε) * p.eta_PQ * ProdPs +
               (up_nh4Pl + nitrate_oxygen_ratio * up_no3Pl + ε) / (up_nPl + ε) * p.eta_PQ * ProdPl -
               f_O2resp * p.eta_O2resp * excrmin -
               (f_O2resp * p.eta_O2resp + f_SO4resp * p.eta_SO4resp) * remtot -
               redfield_N_C * p.eta_nitrif * nitrif) * mol_per_gram_C

    # carbonate dynamics: ΔTA = ΔNH4 - ΔNO3 - ΔPO4 + 2ΔCO3, plus ΔH2S = -ΔO2/2 when O2 ≤ 0
    rhs_dic = mol_per_gram_C * (excrmin + remtot - Prod - rhs_cal)
    rhs_alk = mol_per_gram_C * ((rhs_amm - rhs_nit) * redfield_N_C - rhs_pho * redfield_P_C - 2 * rhs_cal) -
              rhs_oxy * (1 - f_O2pos) / 2

    rates = (no3 = rhs_nit, nh4 = rhs_amm, pho = rhs_pho, sil = rhs_sil, oxy = rhs_oxy,
             fla = rhs_fla, dia = rhs_dia, flachl = rhs_flachl, diachl = rhs_diachl,
             microzoo = rhs_microzoo, mesozoo = rhs_mesozoo, det = rhs_det, opa = rhs_opa,
             cal = rhs_cal, dom = rhs_dom, dic = rhs_dic, alk = rhs_alk)

    diagnostics = (primprod = Prod,
                   secprod = prodZs + prodZl,
                   NlimPs = up_nPs,
                   NlimPl = up_nPl,
                   PlimPs = up_phoPs,
                   PlimPl = up_phoPl,
                   SilimPl = up_silPl,
                   NutlimPs = NutlimPs,
                   NutlimPl = NutlimPl,
                   denit = denit * redfield_N_C * mol_per_gram_C,
                   c2chl_fla = 1 / chl2c_fla,
                   c2chl_dia = 1 / chl2c_dia,
                   cocco = calc_Ps * fla,
                   coccochl = calc_Ps * flachl)

    return (; rates, diagnostics)
end

# 1 at or above `threshold`, 0 below.
@inline at_least(x, threshold) = ifelse(x >= threshold, one(x), zero(x))

# The temperature factor of a rate given at 10 °C.
@inline q10_factor(q10, temp) = q10^((temp - 10) / 10)

"""
    redox_switches(oxy, no3)

The unitary switches `(f_O2pos, f_O2resp, f_denit, f_SO4resp)`: oxic respiration while there is
oxygen, then denitrification while there is nitrate, then sulphate reduction.
"""
@inline function redox_switches(oxy, no3)
    on, off = one(oxy), zero(oxy)
    if oxy <= 0
        return no3 > 0 ? (off, off, on, off) : (off, off, off, on)
    else
        return (on, on, off, off)
    end
end

"""
    oxygen_saturation_concentration(constants, temp, salt)

Oxygen saturation in mmol O₂ m⁻³ at temperature `temp` (°C) and salinity `salt`, after Weiss (1970).
"""
@inline function oxygen_saturation_concentration(constants, temp, salt)
    A1, A2, A3, A4 = constants.oxygen_saturation.temperature
    T = (temp + constants.zero_celsius) / 100
    log_saturation = A1 + A2 / T + A3 * log(T) + A4 * T + salt * evalpoly(T, constants.oxygen_saturation.salinity)
    return exp(log_saturation) * constants.gas_molar_density
end

"""
    transfer_velocity(constants, schmidt, temp, wind)

The Wanninkhof (2014) gas transfer velocity in m s⁻¹ at temperature `temp` (°C, clamped to
-2–40 °C) and 10 m wind `wind` (m s⁻¹), for a gas whose Schmidt number is the polynomial `schmidt`.
"""
@inline function transfer_velocity(constants, schmidt, temp, wind)
    sc = evalpoly(clamp(temp, -2, 40), schmidt)
    return constants.gas_transfer * max(wind, 0)^2 / sqrt(sc / 660)
end

"""
    oxygen_surface_flux(constants, oxy, temp, salt, wind)

The air-sea oxygen flux of `niva_ecosmo`'s `do_surface` in mmol O₂ m⁻² s⁻¹, positive into the
water: the Wanninkhof (2014) transfer velocity times the departure from saturation.
"""
@inline function oxygen_surface_flux(constants, oxy, temp, salt, wind)
    saturation = oxygen_saturation_concentration(constants, temp, max(salt, 0))
    return transfer_velocity(constants, constants.schmidt_O2, temp, wind) * (saturation - oxy)
end

"""
    bottom_fluxes(p::EcosmoParameters, constants, c, s, temp, om_cal, tbs)

`niva_ecosmo`'s `do_bottom` at one water column: the fluxes from the sediment into the bottom
cell, per unit area and positive upward into the water (`fluxes`), and the sources of the four
sediment pools (`sediment`).

`c` holds the bottom-cell concentrations (`oxy`, `no3`, `det`, `opa`, `cal`), `s` the sediment
pools (`sed1` organic matter, `sed2` opal, `sed3` adsorbed phosphate, `sed4` calcite), and `tbs`
is the bottom stress in Pa.
"""
@inline function bottom_fluxes(p::EcosmoParameters, constants, c, s, temp, om_cal, tbs)
    (; oxy, det, opa, cal) = c
    (; sed1, sed2, sed3, sed4) = s
    (; redfield_N_C, redfield_P_C, mol_per_gram_C) = constants

    om_cal = max(om_cal, 0)

    f_O2pos, f_O2resp, f_denit, f_SO4resp = redox_switches(oxy, c.no3)

    # resuspension above the critical bottom stress, settling at the sinking speed below it
    resuspending = tbs >= p.crBotStr
    Rsd = ifelse(resuspending, min(p.resuspRt, p.resuspRt * tbs^2 * 100), zero(p.resuspRt))
    Rds = ifelse(resuspending, zero(p.sinkD), p.sinkD)
    RdsO = ifelse(resuspending, zero(p.sinkOPAL), p.sinkOPAL)
    RdsC = ifelse(resuspending, zero(p.sinkCAL), p.sinkCAL)

    # sediment organic matter (sed1)
    settle = Rds * det
    resusp = Rsd * sed1
    remSED = p.reminSED * exp(p.aTreminSED * temp) * sed1
    burSED = (p.burialRt + p.burialRt2 * sed1) * sed1
    denitSED = f_denit * p.eta_denit * remSED
    lossND = f_O2pos * oxy / (p.KO2nitrifSED + oxy) * p.gammaND * remSED  # coupled nitrification-denitrification

    # sediment opal (sed2)
    settleO = RdsO * opa
    resuspO = Rsd * sed2
    dissSEDO = p.dissSEDO * sed2
    burSEDO = (p.burialRtO + p.burialRt2O * sed2) * sed2

    # sediment phosphate (sed3)
    yt2 = max(oxy / 375, 0)
    yt1 = yt2^2 / (p.RelSEDp2^2 + yt2^2)
    releaseP = p.releaseP * exp(p.aTreleaseP * temp) * (1 - p.RelSEDp1 * yt1) * sed3

    # sediment calcite (sed4)
    settleC = RdsC * cal
    resuspC = Rsd * sed4
    rdissC = max(p.dissSEDCmax * max(1 - om_cal, 0)^p.ndissSEDC, p.dissSEDCmin)
    dissSEDC = rdissC * sed4
    burSEDC = (p.burialRtC + p.burialRt2C * sed4) * sed4

    flux_oxy = -((f_O2resp * p.eta_O2resp + f_SO4resp * p.eta_SO4resp) * remSED +
                 redfield_N_C * p.eta_nitrif * lossND) * mol_per_gram_C

    flux_alk = mol_per_gram_C * ((remSED - lossND + denitSED) * redfield_N_C - releaseP * redfield_P_C + 2 * dissSEDC) -
               flux_oxy * (1 - f_O2pos) / 2

    fluxes = (det = resusp - settle,
              no3 = -denitSED,
              nh4 = remSED - lossND,
              opa = resuspO - settleO,
              sil = dissSEDO,
              pho = releaseP,
              cal = resuspC - settleC,
              oxy = flux_oxy,
              dic = mol_per_gram_C * (remSED + dissSEDC),
              alk = flux_alk)

    sediment = (sed1 = settle - resusp - remSED - burSED,
                sed2 = settleO - resuspO - dissSEDO - burSEDO,
                sed3 = remSED - releaseP,
                sed4 = settleC - resuspC - dissSEDC - burSEDC)

    return (; fluxes, sediment)
end

"""
    sinking_speeds(p::EcosmoParameters, depth)

The sinking speeds of detritus, opal and calcite at `depth` (m, positive down), in m s⁻¹ and
positive downward: the minimum speed plus a linear increase with depth (`get_vertical_movement`).
"""
@inline function sinking_speeds(p::EcosmoParameters, depth)
    depth = abs(depth)
    return (det = p.sinkD + p.cz_sinkD * depth,
            opa = p.sinkOPAL + p.cz_sinkOPAL * depth,
            cal = p.sinkCAL + p.cz_sinkCAL * depth)
end
