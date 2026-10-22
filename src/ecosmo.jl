# The ECOSMO processes as pointwise functions of one cell, ported from `niva_ecosmo.F90`.
#
# These are what both drivers call: the OceanBioME tendency kernels and the box model. Every
# expression keeps the Fortran evaluation order — `a*b/c*d` is `((a*b)/c)*d` in both languages,
# and Julia does not contract into fused multiply-adds unless asked to — so that the two
# implementations round alike. Fortran's integer powers (`x**2`, `x**3`) are written as the
# multiplications gfortran expands them into, and `x**2.0_rk` likewise, since GCC folds it to
# `x*x`.

"""
    interior_rates(p::EcosmoParameters, c, temp, D, E0, om_cal)

The interior sources of every ECOSMO tracer and the diagnostics of `niva_ecosmo`'s `do`, at one point.

`c` holds the concentrations by FABM name (`no3`, `nh4`, `pho`, `sil`, `oxy`, `fla`, `dia`,
`flachl`, `diachl`, `microzoo`, `mesozoo`, `det`, `opa`, `cal`, `dom`), `temp` is in °C, `D` the
day length as a fraction of the day, `E0` the scalar PAR in W m⁻² and `om_cal` the calcite
saturation state.

Returns `(; rates, diagnostics)`, `rates` per second and including `dic` and `alk`, which the
Fortran adds to the carbonate module's state variables when `couple_co2` is on.
"""
@inline function interior_rates(p::EcosmoParameters, c, temp, D, E0, om_cal)
    no3 = c.no3
    nh4 = c.nh4
    pho = c.pho
    sil = c.sil
    oxy = c.oxy
    fla = c.fla
    dia = c.dia
    flachl = c.flachl
    diachl = c.diachl
    microzoo = c.microzoo
    mesozoo = c.mesozoo
    det = c.det
    opa = c.opa
    cal = c.cal
    dom = c.dom

    E0c = max(0.0, E0)

    # Loss terms are switched off below a threshold biomass, to keep plankton from extinction.
    # `max(sign(-1.0_rk, x - 0.1_rk), 0.0_rk)`: 1 when x ≥ 0.1, 0 below.
    fla_loss = max(copysign(1.0, fla - 0.1), 0.0)
    dia_loss = max(copysign(1.0, dia - 0.1), 0.0)
    mic_loss = max(copysign(1.0, microzoo - 0.01), 0.0)
    mes_loss = max(copysign(1.0, mesozoo - 0.01), 0.0)

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
    t_silPl = max(sil - p.minSiPl, 0.0)
    up_silPl = t_silPl / (p.rSiPl + t_silPl)
    NutlimPs = min(up_nPs, up_phoPs)
    NutlimPl = min(up_nPl, up_phoPl, up_silPl)

    # maximum phytoplankton growth rates, corrected for temperature and fractional day length
    muPsc = p.q10Ps^((temp - 10.0) / 10.0) * (p.gammaDPs + 0.5) / (p.gammaDPs + D) * p.muPs
    muPlc = p.q10Pl^((temp - 10.0) / 10.0) * (p.gammaDPl + 0.5) / (p.gammaDPl + D) * p.muPl

    # further correction for light limitation
    ChlCPs = flachl / max(fla, 1e-8)
    ChlCPs = max(p.MINchl2cPs, ChlCPs)
    ChlCPs = min(p.MAXchl2cPs, ChlCPs)
    ChlCPl = diachl / max(dia, 1e-8)
    ChlCPl = max(p.MINchl2cPl, ChlCPl)
    ChlCPl = min(p.MAXchl2cPl, ChlCPl)

    if E0 > 1e-8
        muPsc = muPsc * (1.0 - exp(-p.alfaPs * E0 * ChlCPs / muPsc)) * exp(-p.betaPs * E0 * ChlCPs / muPsc)
        muPlc = muPlc * (1.0 - exp(-p.alfaPl * E0 * ChlCPl / muPlc)) * exp(-p.betaPl * E0 * ChlCPl / muPlc)

        # ratio of chlorophyll-a synthesis to carbon fixation rates
        chl2c_fla = p.MAXchl2cPs * (muPsc * fla) / (p.alfaPs * E0 * flachl)
        chl2c_dia = p.MAXchl2cPl * (muPlc * dia) / (p.alfaPl * E0 * diachl)
        chl2c_fla = max(p.MINchl2cPs, chl2c_fla)
        chl2c_fla = min(p.MAXchl2cPs, chl2c_fla)
        chl2c_dia = max(p.MINchl2cPl, chl2c_dia)
        chl2c_dia = min(p.MAXchl2cPl, chl2c_dia)
    else
        # The E0 → 0 limit of the expressions above (L'Hôpital).
        muPsc = 0.0
        muPlc = 0.0
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
    NutfacmPs = p.amortPs / (1.0 + NutlimPs * (p.amortPs - 1.0))
    NutfacmPl = p.amortPl / (1.0 + NutlimPl * (p.amortPl - 1.0))
    mort1Ps = NutfacmPs * p.m1Ps * fla_loss * fla
    mort1Pl = NutfacmPl * p.m1Pl * dia_loss * dia
    mort2Ps = NutfacmPs * p.m2Ps * fla * fla_loss * fla
    mort2Pl = NutfacmPl * p.m2Pl * dia * dia_loss * dia
    mortPs = mort1Ps + mort2Ps
    mortPl = mort1Pl + mort2Pl

    # grazing/feeding fluxes
    Fs = p.prefZsPs * fla + p.prefZsPl * dia + p.prefZsZs * microzoo + p.prefZsD * det
    Fl = p.prefZlPs * fla + p.prefZlPl * dia + p.prefZlZs * microzoo + p.prefZlZl * mesozoo + p.prefZlD * det

    etZs = p.q10Zs^((temp - 10.0) / 10.0)
    etZl = p.q10Zl^((temp - 10.0) / 10.0)

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
    egesZs = (1.0 - p.gammaZsP) * (ZsonPs + ZsonPl + ZsonZs) + (1.0 - p.gammaZsD) * ZsonD
    excrZs = p.excZs * etZs * mic_loss * microzoo + p.fAexcZs * (p.gammaZsP * (ZsonPs + ZsonPl + ZsonZs) + p.gammaZsD * ZsonD)
    etmZs = p.q10mZs^((temp - 10.0) / 10.0)
    mortZs = p.mZs * etmZs * mic_loss * microzoo
    prodZs = ZsonPs + ZsonPl + ZsonZs + ZsonD - egesZs
    rhs_microzoo = prodZs - excrZs - mortZs - ZsonZs - ZlonZs

    # mesozooplankton, ingestion corrected for sloppy feeding
    inZlPs = (1.0 - p.fGslpZlPs) * ZlonPs
    inZlPl = (1.0 - p.fGslpZlPl) * ZlonPl
    inZlZs = (1.0 - p.fGslpZlZs) * ZlonZs
    inZlZl = (1.0 - p.fGslpZlZl) * ZlonZl
    inZlPZ = inZlPs + inZlPl + inZlZs + inZlZl
    inZlD = (1.0 - p.fGslpZlD) * ZlonD
    egesZl = (1.0 - p.gammaZlP) * inZlPZ + (1.0 - p.gammaZlD) * inZlD
    excrZl = p.excZl * etZl * mes_loss * mesozoo + p.fAexcZl * (p.gammaZlP * inZlPZ + p.gammaZlD * inZlD)
    etmZl = p.q10mZl^((temp - 10.0) / 10.0)
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
    temp² = temp * temp
    frem = p.reminD * (1.0 + 20.0 * (temp² / (169.0 + temp²)))  # 13._rk**2 folds to 169
    remdet = frem * det
    rhs_det = (dxxt - dxxdom) - ZsonD - ZlonD - remdet

    # labile dissolved organic matter
    slpdom = (p.fGslpZlPs * ZlonPs + p.fGslpZlPl * ZlonPl + p.fGslpZlZs * ZlonZs + p.fGslpZlZl * ZlonZl +
              p.fGslpZlD * ZlonD)
    excrt = excrZs + excrZl
    excrdom = p.fexcdomZs * excrZs + p.fexcdomZl * excrZl
    excrmin = excrt - excrdom
    fremDOM = 10.0 * frem
    remdom = fremDOM * dom
    rhs_dom = dxxdom + slpdom + excrdom - remdom

    # nitrate
    nitrif = nh4 * p.nitrifmax * exp(p.aTnitrif * temp) *          # T-dependence
             f_O2pos * oxy / (p.KO2nitrif + oxy) *                  # O2-dependence
             1.0 / (1.0 + p.aE0nitrif * E0c^p.bE0nitrif)            # PAR-dependence
    remtot = remdet + remdom
    denit = f_denit * p.eta_denit * remtot
    rhs_nit = -(up_no3Ps + 0.5e-10) / (up_nPs + 1.0e-10) * ProdPs -
              (up_no3Pl + 0.5e-10) / (up_nPl + 1.0e-10) * ProdPl +
              nitrif -
              denit

    # ammonium
    rhs_amm = -(up_nh4Ps + 0.5e-10) / (up_nPs + 1.0e-10) * ProdPs -
              (up_nh4Pl + 0.5e-10) / (up_nPl + 1.0e-10) * ProdPl +
              excrmin +
              remtot -
              nitrif

    # phosphate
    rhs_pho = -Prod + excrmin + remtot

    # silicate: opal dissolution after Ridgwell et al. (2002)
    sil_eq = 10.0^(6.44 - 968.0 / (temp + 273.15))
    u_opal = min(1.0, max(0.0, (sil_eq - sil * REDFIELD_SI_C * MOL_PER_GRAM_C) / sil_eq))
    ndissO = 0.225 * (1.0 + temp / 15.0) * u_opal + 0.775 * fortran_pow37(1.0 + temp / 400.0) * u_opal^9.25
    dissO = p.dissO * ndissO * opa
    rhs_sil = -ProdPl + dissO

    # opal
    rhs_opa = mortPl + ZsonPl + ZlonPl - dissO

    # calcite: rain ratio within calcifiers, scaled to the small phytoplankton as a whole
    RainR = p.Rain0 * max(0.0, (om_cal - 1.0) / (om_cal - 1.0 + p.Kcalom))
    t = max(0.0, temp)
    calc_Ps = min((1.0 - up_phoPs), up_nPs) * (t / (2.0 + t))  # calcifying fraction of Ps
    RainR = RainR * calc_Ps
    RainR = max(RainR, 0.005)
    calc = RainR * (((1.0 - p.frmort1Ps) * mort1Ps + (1.0 - p.frmort2Ps) * mort2Ps) +
                    ZsonPs * (1.0 - p.fdissCZs) +
                    ZlonPs * (1.0 - p.fdissCZl))
    diss = p.dissCmax * max(1.0 - om_cal, 0.0)^p.ndissC * cal
    rhs_cal = calc - diss

    # oxygen
    rhs_oxy = ((up_nh4Ps + 1.226 * up_no3Ps + 1.0e-10) / (up_nPs + 1.0e-10) * p.eta_PQ * ProdPs +
               (up_nh4Pl + 1.226 * up_no3Pl + 1.0e-10) / (up_nPl + 1.0e-10) * p.eta_PQ * ProdPl -
               f_O2resp * p.eta_O2resp * excrmin -
               (f_O2resp * p.eta_O2resp + f_SO4resp * p.eta_SO4resp) * remtot -
               REDFIELD_N_C * p.eta_nitrif * nitrif) * MOL_PER_GRAM_C

    # carbonate dynamics: ΔTA = ΔNH4 - ΔNO3 - ΔPO4 + 2ΔCO3, plus ΔH2S = -ΔO2/2 when O2 ≤ 0
    rhs_dic = MOL_PER_GRAM_C * (excrmin + remtot - Prod - rhs_cal)
    rhs_alk = MOL_PER_GRAM_C * ((rhs_amm - rhs_nit) * REDFIELD_N_C -
                                rhs_pho * REDFIELD_P_C -
                                2.0 * rhs_cal) -
              0.5 * rhs_oxy * (1.0 - f_O2pos)

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
                   denit = denit * REDFIELD_N_C * MOL_PER_GRAM_C,
                   c2chl_fla = 1.0 / chl2c_fla,
                   c2chl_dia = 1.0 / chl2c_dia,
                   cocco = calc_Ps * fla,
                   coccochl = calc_Ps * flachl)

    return (; rates, diagnostics)
end

"""
    redox_switches(oxy, no3)

The unitary switches `(f_O2pos, f_O2resp, f_denit, f_SO4resp)`: oxic respiration while there is
oxygen, then denitrification while there is nitrate, then sulphate reduction.
"""
@inline function redox_switches(oxy, no3)
    if oxy <= 0.0
        return no3 > 0.0 ? (0.0, 0.0, 1.0, 0.0) : (0.0, 0.0, 0.0, 1.0)
    else
        return (1.0, 1.0, 0.0, 0.0)
    end
end

"""
    fortran_pow37(x)

`x**37` as gfortran computes it: the multiplication chain GCC's `powi` expansion builds from its
`powi_table`, x³⁷ = x³³ · x⁴ with x³³ = x¹⁷ · x¹⁶ and x¹⁷ = x¹⁶ · x. Julia's own `x^37` is more
accurate, which is exactly why it would round differently.
"""
@inline function fortran_pow37(x)
    x² = x * x
    x⁴ = x² * x²
    x⁸ = x⁴ * x⁴
    x¹⁶ = x⁸ * x⁸
    x¹⁷ = x¹⁶ * x
    x³³ = x¹⁷ * x¹⁶
    return x³³ * x⁴
end

"""
    oxygen_saturation_concentration(ETW, X1X)

Oxygen saturation in mmol O₂ m⁻³ at temperature `ETW` (°C) and salinity `X1X`, after Weiss (1970).
"""
@inline function oxygen_saturation_concentration(ETW, X1X)
    A1 = -173.4292
    A2 = 249.6339
    A3 = 143.3483
    A4 = -21.8492
    B1 = -0.033096
    B2 = 0.014259
    B3 = -0.0017
    T = 273.15

    ABT = ETW + T
    ABT100 = ABT / 100.0

    OSAT = A1 + A2 * (100.0 / ABT) + A3 * log(ABT / 100.0) +
           A4 * (ABT / 100.0) +
           X1X * (B1 + B2 * (ABT / 100.0) + B3 * (ABT100 * ABT100))

    OSAT = exp(OSAT)
    return OSAT * 1000.0 / IDEAL_GAS_VOLUME
end

# Volume of an ideal gas at 0 °C and 1 atm, `(R * 273.15_rk / P) * 1000._rk`.
const IDEAL_GAS_VOLUME = (8.3145 * 273.15 / 101325.0) * 1000.0

"""
    oxygen_surface_flux(oxy, temp, salt, wind)

The air-sea oxygen flux of `niva_ecosmo`'s `do_surface` in mmol O₂ m⁻² s⁻¹, positive into the
water: the Wanninkhof (2014) transfer velocity times the departure from saturation.
"""
@inline function oxygen_surface_flux(oxy, temp, salt, wind)
    X1X = max(salt, 0.0)
    wnd = max(wind, 0.0)

    OSAT = oxygen_saturation_concentration(temp, X1X)

    # Schmidt number for O2 following Wanninkhof (2014)
    T = max(min(temp, 40.0), -2.0)
    T² = T * T
    sc = 1920.4 - 135.6 * T + 5.2122 * T² - 0.10939 * (T² * T) + 0.00093777 * (T² * T²)

    ko2o = 0.251 * (wnd * wnd) * (sc / 660.0)^(-0.5)  # cm/hr
    ko2o = ko2o / 360000.0                              # m/s

    return ko2o * (OSAT - oxy)
end

"""
    bottom_fluxes(p::EcosmoParameters, c, s, temp, om_cal, tbs)

`niva_ecosmo`'s `do_bottom` at one water column: the fluxes from the sediment into the bottom
cell, per unit area and positive upward into the water (`fluxes`), and the sources of the four
sediment pools (`sediment`).

`c` holds the bottom-cell concentrations (`oxy`, `no3`, `det`, `opa`, `cal`), `s` the sediment
pools (`sed1` organic matter, `sed2` opal, `sed3` adsorbed phosphate, `sed4` calcite), and `tbs`
is the bottom stress in Pa.
"""
@inline function bottom_fluxes(p::EcosmoParameters, c, s, temp, om_cal, tbs)
    oxy = c.oxy
    det = c.det
    opa = c.opa
    cal = c.cal
    sed1 = s.sed1
    sed2 = s.sed2
    sed3 = s.sed3
    sed4 = s.sed4

    om_cal = max(om_cal, 0.0)

    f_O2pos, f_O2resp, f_denit, f_SO4resp = redox_switches(oxy, c.no3)

    # resuspension above the critical bottom stress, settling at the sinking speed below it
    if tbs >= p.crBotStr
        Rsd = min(p.resuspRt, p.resuspRt * (tbs * tbs) * 100.0)
        Rds = 0.0
        RdsO = 0.0
        RdsC = 0.0
    else
        Rsd = 0.0
        Rds = p.sinkD
        RdsO = p.sinkOPAL
        RdsC = p.sinkCAL
    end

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
    yt2 = max(oxy / 375.0, 0.0)
    yt1 = (yt2 * yt2) / (p.RelSEDp2 * p.RelSEDp2 + yt2 * yt2)
    releaseP = p.releaseP * exp(p.aTreleaseP * temp) * (1.0 - p.RelSEDp1 * yt1) * sed3

    # sediment calcite (sed4)
    settleC = RdsC * cal
    resuspC = Rsd * sed4
    rdissC = p.dissSEDCmax * (max(1.0 - om_cal, 0.0))^p.ndissSEDC
    rdissC = max(rdissC, p.dissSEDCmin)
    dissSEDC = rdissC * sed4
    burSEDC = (p.burialRtC + p.burialRt2C * sed4) * sed4

    flux_oxy = -((f_O2resp * p.eta_O2resp + f_SO4resp * p.eta_SO4resp) * remSED +
                 REDFIELD_N_C * p.eta_nitrif * lossND) * MOL_PER_GRAM_C

    flux_alk = MOL_PER_GRAM_C * ((remSED - lossND + denitSED) * REDFIELD_N_C -
                                 releaseP * REDFIELD_P_C +
                                 2.0 * dissSEDC) -
               0.5 * flux_oxy * (1.0 - f_O2pos)

    fluxes = (det = resusp - settle,
              no3 = -denitSED,
              nh4 = remSED - lossND,
              opa = resuspO - settleO,
              sil = dissSEDO,
              pho = releaseP,
              cal = resuspC - settleC,
              oxy = flux_oxy,
              dic = MOL_PER_GRAM_C * (remSED + dissSEDC),
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
