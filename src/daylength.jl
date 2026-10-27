# `niva_ecosmo_daylength.F90`: the day length from the NOAA declination formula.

"""
    day_length(constants, latitude, day_of_year)

The day length as a fraction of the day at `latitude` (degrees north), on the day containing
`day_of_year` — FABM's `number_of_days_since_start_of_the_year`, zero at midnight on 1 January.

Leap years are not considered, as in the Fortran.
"""
@inline function day_length(constants, latitude, day_of_year)
    θ = constants.radians_per_day * (unsafe_trunc(Int, day_of_year) + 1)
    d0, d1, d2, d3, d4, d5, d6 = constants.declination
    declination = d0 + d1 * cos(θ) + d2 * sin(θ) + d3 * cos(2θ) + d4 * sin(2θ) + d5 * cos(3θ) + d6 * sin(3θ)

    φ = deg2rad(latitude)
    cos_hour_angle = constants.cos_sunrise_zenith / (cos(φ) * cos(declination)) - tan(φ) * tan(declination)

    return acos(clamp(cos_hour_angle, -1, 1)) / π
end
