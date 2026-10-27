# The fixed numbers of the model, as opposed to the `fabm.yaml` parameters of `parameters.jl`.
#
# A process function computes with parameters, exact integers and the entries of `ECOSMO_CONSTANTS`
# only. The table is written once, in `Float64`; `ecosmo_constants(FT)` converts it to the model's
# float type when the model is built, so the arithmetic is in `FT` throughout.

# Redfield ratios and the molar mass of carbon (`redf` in `niva_ecosmo.F90`), which also convert
# units when a model is set up.
const REDFIELD_C_N = 6.625
const REDFIELD_C_P = 106.0
const REDFIELD_C_SI = 6.625
const CARBON_MOLAR_MASS = 12.01  # g C per mol C

const ECOSMO_CONSTANTS = (
    # stoichiometry
    redfield_N_C = 1 / REDFIELD_C_N,
    redfield_P_C = 1 / REDFIELD_C_P,
    redfield_Si_C = 1 / REDFIELD_C_SI,
    mol_per_gram_C = 1 / CARBON_MOLAR_MASS,
    nitrate_oxygen_ratio = 1.226,  # O₂ released by growth on nitrate, relative to growth on ammonium

    # plankton
    phytoplankton_loss_threshold = 0.1,   # mgC m⁻³: no mortality or grazing below
    zooplankton_loss_threshold = 0.01,    # mgC m⁻³
    reference_day_length = 0.5,           # fraction of the day at which day length does not affect growth
    minimum_light = 1e-8,                 # W m⁻²: no photosynthesis below
    minimum_carbon = 1e-8,                # mgC m⁻³: floor of the carbon in the Chl:C ratio
    uptake_regularization = 1e-10,        # keeps the NO₃/NH₄ uptake split defined when both are zero
    minimum_rain_ratio = 0.005,           # calcite to small-phytoplankton carbon

    # opal dissolution (Ridgwell et al. 2002)
    opal_solubility = (6.44, 968.0),  # log10 of the equilibrium silicate: a - b / T(K)
    opal_dissolution = (linear = 0.225, power = 0.775, exponent = 9.25),

    # oxygen saturation (Weiss 1970), ml l⁻¹ from T/100 in K
    oxygen_saturation = (temperature = (-173.4292, 249.6339, 143.3483, -21.8492),
                         salinity = (-0.033096, 0.014259, -0.0017)),
    gas_molar_density = 101325 / (8.3145 * 273.15),  # mol m⁻³ at 0 °C and 1 atm: ml l⁻¹ → mmol m⁻³

    # air-sea gas transfer (Wanninkhof 2014): Schmidt numbers as polynomials in °C, and the
    # transfer velocity 0.251 cm h⁻¹ per (m s⁻¹)² in m s⁻¹
    schmidt_O2 = (1920.4, -135.6, 5.2122, -0.10939, 0.00093777),
    schmidt_CO2 = (2116.8, -136.25, 4.7353, -0.092307, 0.0007555),
    gas_transfer = 0.251 / 360000,

    zero_celsius = 273.15,  # K

    # carbonate system (`niva_ecosmo_carbonate.F90`, the ERSEM carbonate module)
    gas_constant = 83.131,                # cm³ bar mol⁻¹ K⁻¹
    boron_per_salinity = 0.0004128 / 35,  # mol kg⁻¹
    calcium = 0.01028,                    # mol kg⁻¹
    ionic_strength = (19.924, 1.005),     # I = a S / (1000 - b S)
    chlorinity_per_salinity = 1 / 1.80655,
    sulfate_per_chlorinity = 0.14 / 96.065,
    fluoride_per_chlorinity = 0.000067 / 19.9984,
    # ln ks, bisulfate (Dickson 1990): each triple is a / T + b + c ln T. The Fortran writes 771.54
    # without a kind suffix, so it is a single-precision number there, and `771.54f0` here.
    bisulfate = (pure = (-4276.1, 141.328, -23.093),
                 root_ionic = (-13856.0, 324.57, -47.986),
                 ionic = (35474.0, -771.54f0, 114.723),
                 ionic_1p5 = -2698.0,  # times I^1.5 / T
                 ionic_2 = 1776.0,     # times I² / T
                 salinity = 0.001005), # ln(1 - c S)
    # ln kf, hydrogen fluoride (Perez and Fraga 1987): a / T + b + c √S
    fluoride = (874.0, -9.68, 0.111),
    # ln k0, CO₂ solubility (Weiss 1974), T/100 in K: a / T + b + c ln T + S (d + e T + f T²)
    co2_solubility = (temperature = (93.4517, -60.2409, 23.3585),
                      salinity = (0.023517, -0.023656, 0.0047036)),
    co2_partial_molar_volume = 32.3,  # cm³ mol⁻¹
    standard_pressure = 1.01325,      # bar
    # ln kb, boric acid (Dickson 1990), each a polynomial in √S: A / T + B + C ln T + D √S T
    boric_acid = (inverse_T = (-8966.9, -2890.53, -77.942, 1.728, -0.0996),
                  constant = (148.0248, 137.1942, 1.62142),
                  log_T = (-24.4344, -25.085, -0.2474),
                  T = 0.053105),
    # pK1 and pK2 of carbonic acid for pH scales -1, 0 and 1 (Millero 1995 old seawater scale,
    # Millero 2010 seawater and total scales), as (a0, a1, a2, b1, b2, b3, c1, c2, d) in
    # a0 + a1/T + a2 ln T + b1 √S + b2 S + b3 S² + (c1 √S + c2 S)/T + d √S ln T. 5.7789 is
    # single precision in the Fortran, like 771.54 above.
    carbonic_acid_1 = ((-62.008, 3670.7, 9.7944, 0.0, -0.0118, 0.000116, 0.0, 0.0, 0.0),
                       (-126.34048, 6320.813, 19.568224, 13.4038, 0.03206, -0.00005242, -530.659, -5.8210, -2.0664),
                       (-126.34048, 6320.813, 19.568224, 13.4051, 0.03185, -0.00005218, -531.095, -5.7789f0, -2.0663)),
    carbonic_acid_2 = ((4.777, 1394.7, 0.0, 0.0, -0.0184, 0.000118, 0.0, 0.0, 0.0),
                       (-90.18333, 5143.692, 14.613358, 21.3728, 0.1218, -0.0003688, -788.289, -19.189, -3.374),
                       (-90.18333, 5143.692, 14.613358, 21.5724, 0.1212, -0.0003714, -798.292, -18.951, -3.403)),
    # Pressure corrections (Millero 1995): the partial molar volume (cm³ mol⁻¹) and compressibility
    # (cm³ mol⁻¹ bar⁻¹) changes as polynomials in °C. As in the Fortran, the constant compressibility
    # terms of bisulfate and fluoride are not divided by 1000, unlike the others.
    bisulfate_pressure = (volume = (-18.03, 0.0466, 0.000316), compressibility = (-4.53, 0.00009)),
    fluoride_pressure = (volume = (-9.78, -0.009, -0.000942), compressibility = (-3.91, 0.000054)),
    carbonic_acid_1_pressure = (volume = (-25.5, 0.1271), compressibility = (-3.08e-3, 0.0877e-3)),
    carbonic_acid_2_pressure = (volume = (-15.82, -0.0219), compressibility = (1.13e-3, -0.1475e-3)),
    boric_acid_pressure = (volume = (-29.48, 0.1622, -0.002608), compressibility = (-2.84e-3,)),
    calcite_pressure = (volume = (-48.76, 0.5304), compressibility = (-11.76e-3, 0.3692e-3)),
    aragonite_pressure = (volume = (-46.0, 0.530), compressibility = (-11.76e-3, 0.3692e-3)),
    # log10 Ksp (Mucci 1983): a0 + a1 T + a2 / T + a3 log10 T + (b0 + b1 T + b2 / T) √S + c0 S + c1 S^1.5
    calcite_solubility = (temperature = (-171.9065, -0.077993, 2839.319, 71.595),
                          root_salinity = (-0.77712, 0.0028426, 178.34),
                          salinity = (-0.07711, 0.0041249)),
    aragonite_solubility = (temperature = (-171.945, -0.077993, 2903.293, 71.595),
                            root_salinity = (-0.068393, 0.0017276, 88.135),
                            salinity = (-0.10018, 0.0059415)),
    # HALTAFALL's iteration (`CO2CLC`)
    initial_hydrogen_ion = 1e-8,  # mol kg⁻¹
    iteration_tolerance = 1e-5,   # relative to DIC, alkalinity and boron
    secant_threshold = 0.4,       # step in ln[H⁺] below which the iteration turns to the secant

    # day length (`niva_ecosmo_daylength.F90`): the declination's Fourier series in the day angle
    # (Spencer 1971), and the solar zenith angle at sunrise, 90.833°
    declination = (0.006918, -0.399912, 0.070257, -0.006758, 0.000907, -0.002697, 0.001480),
    radians_per_day = 2π / 365,
    cos_sunrise_zenith = cosd(90.833),
)

"""
    ecosmo_constants(FT = Float64)

`ECOSMO_CONSTANTS` in float type `FT`: the fixed numbers the process functions take with the
parameters. `ECOSMO`, `ecosmo_sediment` and `ecosmo_surface_fluxes` build it from their `FT`.
"""
ecosmo_constants(FT::DataType = Float64) = to_float_type(FT, ECOSMO_CONSTANTS)

to_float_type(FT, x::Number) = convert(FT, x)
to_float_type(FT, x::Union{Tuple, NamedTuple}) = map(y -> to_float_type(FT, y), x)
