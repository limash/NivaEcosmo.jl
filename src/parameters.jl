# The ECOSMO parameters, with the names, defaults and units of `niva_ecosmo.F90`'s `initialize`.
#
# FABM converts a parameter once, when it is read: `get_parameter(..., scale_factor = s)` stores
# `value * s`. Rates are given per day and stored per second, and the nutrient half-saturations are
# given in mmol N, P or Si m⁻³ and stored in the mgC m⁻³ the state variables are kept in. The table
# below carries the same scale factor next to each default, so `EcosmoParameters` holds exactly the
# numbers the Fortran model computes with, and a `fabm.yaml` value means the same thing in both.

const SECONDS_PER_DAY = 86400.0                      # `sedy0`
const PER_DAY = 1.0 / SECONDS_PER_DAY                # scale_factor = 1.0_rk/sedy0
const MGC_PER_MMOL_N = REDFIELD_C_N * CARBON_MOLAR_MASS   # scale_factor = redf(1)*redf(6)
const MGC_PER_MMOL_P = REDFIELD_C_P * CARBON_MOLAR_MASS   # scale_factor = redf(2)*redf(6)
const MGC_PER_MMOL_SI = REDFIELD_C_SI * CARBON_MOLAR_MASS # scale_factor = redf(3)*redf(6)

# (name, default, scale factor, units as given, description), in `initialize` order.
const ECOSMO_PARAMETERS = (
    # light attenuation and suspended matter
    (:Exphy, 0.04, 1.0, "m2/mgCHL", "phyto extinction"),
    (:ExPs, 0.0, 1.0, "m2/mgCHL", "small phyto extinction"),
    (:ExPl, 0.0, 1.0, "m2/mgCHL", "large phyto extinction"),
    (:rtsom_cnp, 1.5547, 1.0, "-", "ratio of total suspended organic mass (dry) to total CNP mass"),
    (:rtsim_s, 2.4596, 1.0, "-", "ratio of total suspended inorganic mass (dry) to total silicon mass"),
    (:rtsim_c_mesozoo, 0.1583, 1.0, "-", "ratio of suspended inorganic mass to carbon mass in mesozooplankton"),
    # small phytoplankton
    (:rNH4Ps, 0.20, MGC_PER_MMOL_N, "mmolN/m3", "NH4 half saturation for Ps"),
    (:rNO3Ps, 0.50, MGC_PER_MMOL_N, "mmolN/m3", "NO3 half saturation for Ps"),
    (:rPO4Ps, 0.05, MGC_PER_MMOL_P, "mmolP/m3", "PO4 half saturation for Ps"),
    (:muPs, 1.10, PER_DAY, "1/day", "max growth rate for Ps"),
    (:q10Ps, 1.0, 1.0, "-", "Q_10 temperature coefficient for Ps"),
    (:gammaDPs, 1.0e8, 1.0, "-", "day length adaptation coefficient for Ps"),
    (:alfaPs, 3.127, PER_DAY, "mgC/(mgChl day W/m2)", "initial slope P-I curve for Ps"),
    (:betaPs, 0.0, PER_DAY, "mgC/(mgChl day W/m2)", "photoinhibition parameter for Ps"),
    (:MINchl2cPs, 0.0063, 1.0, "mgChl/mgC", "minimum Chl to C ratio Ps"),
    (:MAXchl2cPs, 0.0481, 1.0, "mgChl/mgC", "maximum Chl to C ratio Ps"),
    (:amortPs, 1.0, 1.0, "-", "ratio max:min Ps mortality under nutrient limitation"),
    (:m1Ps, 0.08, PER_DAY, "1/day", "linear Ps mortality rate"),
    (:m2Ps, 0.0, PER_DAY, "m3/(mgC day)", "quadratic Ps mortality rate"),
    (:frmort1Ps, 0.4, 1.0, "-", "dissolved fraction of linear Ps mortality"),
    (:frmort2Ps, 0.4, 1.0, "-", "dissolved fraction of quadratic Ps mortality"),
    # large phytoplankton
    (:rNH4Pl, 0.20, MGC_PER_MMOL_N, "mmolN/m3", "NH4 half saturation for Pl"),
    (:rNO3Pl, 0.50, MGC_PER_MMOL_N, "mmolN/m3", "NO3 half saturation for Pl"),
    (:rPO4Pl, 0.05, MGC_PER_MMOL_P, "mmolP/m3", "PO4 half saturation for Pl"),
    (:rSiPl, 0.50, MGC_PER_MMOL_SI, "mmolSi/m3", "SiO2 half saturation for Pl"),
    (:minSiPl, 1.0, MGC_PER_MMOL_SI, "mmolSi/m3", "stop Si uptake below this concentration"),
    (:muPl, 1.30, PER_DAY, "1/day", "max growth rate for Pl"),
    (:q10Pl, 1.0, 1.0, "-", "Q_10 temperature coefficient for Pl"),
    (:gammaDPl, 0.5, 1.0, "-", "day length adaptation coefficient for Pl"),
    (:alfaPl, 4.225, PER_DAY, "mgC/(mgChl day W/m2)", "initial slope P-I curve for Pl"),
    (:betaPl, 0.0, PER_DAY, "mgC/(mgChl day W/m2)", "photoinhibition parameter for Pl"),
    (:MINchl2cPl, 0.0063, 1.0, "mgChl/mgC", "minimum Chl to C ratio Pl"),
    (:MAXchl2cPl, 0.0370, 1.0, "mgChl/mgC", "maximum Chl to C ratio Pl"),
    (:amortPl, 1.0, 1.0, "-", "ratio max:min Pl mortality under nutrient limitation"),
    (:m1Pl, 0.04, PER_DAY, "1/day", "linear Pl mortality rate"),
    (:m2Pl, 0.0, PER_DAY, "m3/(mgC day)", "quadratic Pl mortality rate"),
    (:frmort1Pl, 0.4, 1.0, "-", "dissolved fraction of linear Pl mortality"),
    (:frmort2Pl, 0.0, 1.0, "-", "dissolved fraction of quadratic Pl mortality"),
    # small zooplankton
    (:GrZsPs, 1.00, PER_DAY, "1/day", "Grazing rate Zs on Ps"),
    (:GrZsPl, 1.00, PER_DAY, "1/day", "Grazing rate Zs on Pl"),
    (:GrZsZs, 0.00, PER_DAY, "1/day", "Grazing rate Zs on Zs"),
    (:GrZsD, 1.00, PER_DAY, "1/day", "Grazing rate Zs on Det"),
    (:prefZsPs, 0.70, 1.0, "-", "Grazing preference Zs on Ps"),
    (:prefZsPl, 0.25, 1.0, "-", "Grazing preference Zs on Pl"),
    (:prefZsZs, 0.00, 1.0, "-", "Grazing preference Zs on Zs"),
    (:prefZsD, 0.00, 1.0, "-", "Grazing preference Zs on Det"),
    (:RgZs, 40.0, 1.0, "mgC/m3", "Zs grazing half saturation"),
    (:q10Zs, 1.0, 1.0, "-", "Q_10 temperature coefficient for Zs grazing"),
    (:gammaZsP, 0.75, 1.0, "-", "Zs assim. eff. on plankton"),
    (:gammaZsD, 0.75, 1.0, "-", "Zs assim. eff. on det"),
    (:fregesZs, 0.4, 1.0, "-", "dissolved fraction of Zs egestion"),
    (:excZs, 0.08, PER_DAY, "1/day", "Zs basal excretion/respiration rate"),
    (:fAexcZs, 0.0, 1.0, "-", "Zs fraction of assimilated food that is excreted/respired"),
    (:fexcdomZs, 0.0, 1.0, "-", "Zs fraction of excreted/respired matter that is dom"),
    (:mZs, 0.20, PER_DAY, "1/day", "Zs mortality rate"),
    (:q10mZs, 1.0, 1.0, "-", "Q_10 temperature coefficient for Zs mortality"),
    (:frmortZs, 0.4, 1.0, "-", "dissolved fraction of Zs mortality"),
    # large zooplankton
    (:GrZlPs, 0.80, PER_DAY, "1/day", "Grazing rate Zl on Ps"),
    (:GrZlPl, 0.80, PER_DAY, "1/day", "Grazing rate Zl on Pl"),
    (:GrZlZs, 0.50, PER_DAY, "1/day", "Grazing rate Zl on Zs"),
    (:GrZlZl, 0.00, PER_DAY, "1/day", "Grazing rate Zl on Zl"),
    (:GrZlD, 0.80, PER_DAY, "1/day", "Grazing rate Zl on Det"),
    (:prefZlPs, 0.10, 1.0, "-", "Grazing preference Zl on Ps"),
    (:prefZlPl, 0.85, 1.0, "-", "Grazing preference Zl on Pl"),
    (:prefZlZs, 0.15, 1.0, "-", "Grazing preference Zl on Zs"),
    (:prefZlZl, 0.00, 1.0, "-", "Grazing preference Zl on Zl"),
    (:prefZlD, 0.00, 1.0, "-", "Grazing preference Zl on Det"),
    (:RgZl, 40.0, 1.0, "mgC/m3", "Zl grazing half saturation"),
    (:q10Zl, 1.0, 1.0, "-", "Q_10 temperature coefficient for Zl grazing"),
    (:fGslpZlPs, 0.0, 1.0, "-", "Zl fraction of Ps-grazing that is slopped (not ingested)"),
    (:fGslpZlPl, 0.0, 1.0, "-", "Zl fraction of Pl-grazing that is slopped (not ingested)"),
    (:fGslpZlZs, 0.0, 1.0, "-", "Zl fraction of Zs-feeding that is slopped (not ingested)"),
    (:fGslpZlZl, 0.0, 1.0, "-", "Zl fraction of Zl-feeding that is slopped (not ingested)"),
    (:fGslpZlD, 0.0, 1.0, "-", "Zl fraction of Det-feeding that is slopped (not ingested)"),
    (:gammaZlP, 0.75, 1.0, "-", "Zl assim. eff. on plankton"),
    (:gammaZlD, 0.75, 1.0, "-", "Zl assim. eff. on det"),
    (:fregesZl, 0.4, 1.0, "-", "dissolved fraction of Zl egestion"),
    (:excZl, 0.06, PER_DAY, "1/day", "Zl basal excretion/respiration rate"),
    (:fAexcZl, 0.0, 1.0, "-", "Zl fraction of assimilated food that is excreted/respired"),
    (:fexcdomZl, 0.0, 1.0, "-", "Zl fraction of excreted/respired matter that is dom"),
    (:mZl, 0.10, PER_DAY, "1/day", "Zl basal mortality rate"),
    (:mvZl, 0.0, PER_DAY, "1/day", "Zl mortality due to visual predation"),
    (:KE0mvZl, 2.2, 1.0, "W/m2", "Half-saturation PAR for Zl visual predation mortality"),
    (:q10mZl, 1.0, 1.0, "-", "Q_10 temperature coefficient for Zl mortality"),
    (:frmortZl, 0.4, 1.0, "-", "dissolved fraction of Zl mortality"),
    # stoichiometry
    (:eta_PQ, 1.0, 1.0, "-", "photosynthetic quotient: moles O2 produced per mole C fixed"),
    (:eta_O2resp, 1.0, 1.0, "-", "O2 respiration quotient: moles O2 consumed per mole C respired"),
    (:eta_denit, 5.3, 1.0, "-", "denitrification quotient: moles NO3 used per mole N in organic matter"),
    (:eta_SO4resp, 1.0, 1.0, "-", "SO4 respiration quotient: moles negative O2 produced per mole C respired"),
    (:eta_nitrif, 2.0, 1.0, "-", "nitrification quotient: moles O2 consumed per mole NH4 nitrified"),
    # remineralization and regeneration
    (:nitrifmax, 0.1, PER_DAY, "1/day", "maximum nitrification rate"),
    (:aTnitrif, 0.11, 1.0, "1/degC", "temp. control nitrification"),
    (:KO2nitrif, 0.45, 1.0, "mmolO2/m3", "O2 half saturation for nitrification"),
    (:aE0nitrif, 0.0, 1.0, "(W/m2)^-bE0nitrif", "PAR control nitrification"),
    (:bE0nitrif, 1.0, 1.0, "-", "PAR exponent nitrification"),
    (:reminD, 0.003, PER_DAY, "1/day", "detritus remin. rate"),
    (:dissO, 0.015, PER_DAY, "1/day", "opal dissolution rate (0 degC, sil=0)"),
    # calcite formation and dissolution
    (:Rain0, 1.0, 1.0, "-", "maximum Rain Ratio (PIC:POC) within calcifiers"),
    (:Kcalom, 1.0, 1.0, "-", "half-saturation constant for calcifier Rain Ratio dependence on Omega(calcite)"),
    (:fdissCZs, 0.50, 1.0, "-", "fraction of small zooplankton prey calcite that dissolves after ingestion"),
    (:fdissCZl, 0.25, 1.0, "-", "fraction of large zooplankton prey calcite that dissolves after ingestion"),
    (:dissCmax, 0.03, PER_DAY, "1/day", "maximum specific dissolution rate"),
    (:ndissC, 2.22, 1.0, "-", "power of the dissolution law (Keir 1980)"),
    # sinking
    (:sinkD, 5.0, PER_DAY, "m/day", "minimum detritus sinking rate"),
    (:sinkOPAL, 5.0, PER_DAY, "m/day", "minimum OPAL sinking rate"),
    (:sinkCAL, 5.0, PER_DAY, "m/day", "minimum calcite sinking rate"),
    (:cz_sinkD, 0.0, PER_DAY, "1/day", "detritus sinking rate increase per metre below surface"),
    (:cz_sinkOPAL, 0.0, PER_DAY, "1/day", "opal sinking rate increase per metre below surface"),
    (:cz_sinkCAL, 0.0, PER_DAY, "1/day", "calcite sinking rate increase per metre below surface"),
    # sediment
    (:crBotStr, 0.007, 1.0, "N/m2", "critic. bot. stress for resusp."),
    (:resuspRt, 25.0, PER_DAY, "1/day", "resuspension rate"),
    (:burialRt, 1e-5, PER_DAY, "1/day", "linear detritus burial rate"),
    (:burialRt2, 0.0, PER_DAY, "m2/(mgC day)", "quadratic detritus burial parameter"),
    (:reminSED, 0.002, PER_DAY, "1/day", "remineralization rate"),
    (:aTreminSED, 0.15, 1.0, "1/degC", "temp. control remineralization"),
    (:gammaND, 0.5, 1.0, "-", "fraction of oxic remineralization lost to coupled nitrification-denitrification"),
    (:KO2nitrifSED, 4.5, 1.0, "mmolO2/m3", "O2 half saturation for nitrification in sediments"),
    (:dissSEDO, 0.0002, PER_DAY, "1/day", "sed. opal dissolution rate"),
    (:burialRtO, 1e-5, PER_DAY, "1/day", "linear opal burial rate"),
    (:burialRt2O, 0.0, PER_DAY, "m2/(mgC day)", "quadratic opal burial parameter"),
    (:dissSEDCmax, 30.0, PER_DAY, "1/day", "sed. max. dissolution rate calcite"),
    (:dissSEDCmin, 0.05, PER_DAY, "1/day", "sed. min. dissolution rate calcite"),
    (:ndissSEDC, 2.22, 1.0, "-", "sed. power of calcite dissolution law (Keir 1980)"),
    (:burialRtC, 1e-5, PER_DAY, "1/day", "linear calcite burial rate"),
    (:burialRt2C, 0.0, PER_DAY, "m2/(mgC day)", "quadratic calcite burial parameter"),
    (:releaseP, 0.002, PER_DAY, "1/day", "P sed. release rate"),
    (:aTreleaseP, 0.15, 1.0, "1/degC", "temp. control P sed. release"),
    (:RelSEDp1, 0.15, 1.0, "-", "P sed. release O2-dependence max. fractional reduction"),
    (:RelSEDp2, 0.10, 1.0, "-", "P sed. release O2-dependence half-saturation coefficient"),
)

const ECOSMO_PARAMETER_NAMES = map(first, ECOSMO_PARAMETERS)

@eval struct EcosmoParameters{FT}
    $([:($name::FT) for name in ECOSMO_PARAMETER_NAMES]...)
end

"""
    EcosmoParameters{FT}

Every `niva_ecosmo` parameter, in the units the model computes with: rates per second and nutrient
half-saturations in mgC m⁻³. Build one with `EcosmoParameters(; kwargs...)`, which takes the
`fabm.yaml` units, or read a `fabm.yaml` with `fabm_configuration`.

Generated from `ECOSMO_PARAMETERS`, so the field list, the defaults and the unit conversion are
written down once.
"""
EcosmoParameters

"""
    EcosmoParameters(FT = Float64; kwargs...)

The ECOSMO parameters, each keyword in the units `fabm.yaml` gives it (per day, mmol N m⁻³, ...)
and defaulting to `niva_ecosmo.F90`'s own default. Each value is multiplied by its FABM scale factor
on the way in, exactly as `get_parameter` does.

An unknown keyword is an error, because a misspelt parameter would otherwise keep its default with
nothing to say so. FABM rejects one too.
"""
function EcosmoParameters(FT::DataType = Float64; kwargs...)
    unknown = setdiff(keys(kwargs), ECOSMO_PARAMETER_NAMES)
    isempty(unknown) || throw(ArgumentError("unknown ECOSMO parameter(s): $(join(unknown, ", "))"))

    values = map(ECOSMO_PARAMETERS) do (name, default, scale_factor)
        convert(FT, Float64(get(kwargs, name, default)) * scale_factor)
    end

    return EcosmoParameters{FT}(values...)
end

Base.summary(::EcosmoParameters{FT}) where {FT} = "EcosmoParameters{$FT}"
