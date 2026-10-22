# `niva_ecosmo_daylength.F90`: the day length from the NOAA declination formula.

# `cos(90.833_rk*pi/180._rk)` is a constant expression, which gfortran folds at compile time with
# correctly rounded arithmetic, so it is evaluated in extended precision here.
const COS_SUNRISE_ZENITH = Float64(cos(big(90.833 * PI_RK / 180.0)))

"""
    day_length(latitude, day_of_year)

The day length as a fraction of the day at `latitude` (degrees north), on the day containing
`day_of_year` — FABM's `number_of_days_since_start_of_the_year`, zero at midnight on 1 January.

Leap years are not considered, as in the Fortran.
"""
@inline function day_length(latitude, day_of_year)
    iday = unsafe_trunc(Int, day_of_year) + 1

    th0 = PI_RK * iday / 182.5
    th02 = 2.0 * th0
    th03 = 3.0 * th0

    sundec = 0.006918 - 0.399912 * cos(th0) + 0.070257 * sin(th0) -
             0.006758 * cos(th02) + 0.000907 * sin(th02) -
             0.002697 * cos(th03) + 0.001480 * sin(th03)

    latr = latitude * PI_RK / 180.0

    return acos(min(1.0, max(-1.0, COS_SUNRISE_ZENITH / (cos(latr) * cos(sundec)) - tan(latr) * tan(sundec)))) / PI_RK
end
