# Fixed constants of `niva_ecosmo.F90` (`redf`, set in `initialize`) and the other Fortran modules.
#
# Each is written the way the Fortran computes it, so the rounding is the same: `redf(11:20)` are
# `1._rk/redf(1:10)`, not their decimal values.

const REDFIELD_C_N = 6.625          # redf(1),  C:N
const REDFIELD_C_P = 106.0          # redf(2),  C:P
const REDFIELD_C_SI = 6.625         # redf(3),  C:Si
const CARBON_MOLAR_MASS = 12.01     # redf(6),  g C per mol C

const REDFIELD_N_C = 1.0 / REDFIELD_C_N     # redf(11)
const REDFIELD_P_C = 1.0 / REDFIELD_C_P     # redf(12)
const REDFIELD_SI_C = 1.0 / REDFIELD_C_SI   # redf(13)
const MOL_PER_GRAM_C = 1.0 / CARBON_MOLAR_MASS  # redf(16)

const PI_RK = 3.141592653589793     # the `pi` parameter of niva_ecosmo_daylength.F90
