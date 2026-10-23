# Reading a niva-ecosmo `fabm.yaml`, so the Julia model runs with the configuration the Fortran
# one does. Instances are found by their `model`, not by their name, since names vary between
# configurations. FABM matches parameter names without regard to case, and so does this.

const ECOSMO_MODEL = "niva/niva_ecosmo"
const CARBONATE_MODEL = "niva/niva_ecosmo_carbonate"
const LIGHT_MODEL = "gotm/light"

# The Fortran `initial_value`s, written as the Fortran computes them (`5.0_rk*redf(1)*redf(6)` is
# `(5.0*6.625)*12.01`).
const ECOSMO_DEFAULT_INITIAL_VALUES = (
    no3 = 5.0 * REDFIELD_C_N * CARBON_MOLAR_MASS,
    nh4 = 0.1 * REDFIELD_C_N * CARBON_MOLAR_MASS,
    pho = 0.3 * REDFIELD_C_P * CARBON_MOLAR_MASS,
    sil = 5.0 * REDFIELD_C_SI * CARBON_MOLAR_MASS,
    oxy = 85.0,
    fla = 1e-4 * REDFIELD_C_N * CARBON_MOLAR_MASS,
    dia = 1e-4 * REDFIELD_C_N * CARBON_MOLAR_MASS,
    flachl = 1e-4 * REDFIELD_C_N * CARBON_MOLAR_MASS / 20.0,
    diachl = 1e-4 * REDFIELD_C_N * CARBON_MOLAR_MASS / 27.0,
    microzoo = 1e-6 * REDFIELD_C_N * CARBON_MOLAR_MASS,
    mesozoo = 1e-6 * REDFIELD_C_N * CARBON_MOLAR_MASS,
    det = 2.0 * REDFIELD_C_N * CARBON_MOLAR_MASS,
    opa = 2.0 * REDFIELD_C_SI * CARBON_MOLAR_MASS,
    cal = 0.05,
    dom = 3.0 * REDFIELD_C_N * CARBON_MOLAR_MASS,
    dic = 2200.0,
    alk = 2300.0,
    sed1 = 20.0 * REDFIELD_C_N * CARBON_MOLAR_MASS * (1.0 / 14.007),
    sed2 = 20.0 * REDFIELD_C_SI * CARBON_MOLAR_MASS * (1.0 / 28.09),
    sed3 = 2.0 * REDFIELD_C_P * CARBON_MOLAR_MASS * (1.0 / 30.97),
    sed4 = 0.05,
)

"""
    fabm_configuration(path, FT = Float64)

The ECOSMO configuration in the `fabm.yaml` at `path`, as
`(; parameters, pH_scale, light, initial_values)`:

- `parameters`: the `niva/niva_ecosmo` instance's [`EcosmoParameters`](@ref), in `FT`, which must
  be the float type of the grid they are used on;
- `pH_scale`: the `niva/niva_ecosmo_carbonate` instance's `pHscale`;
- `light`: the `gotm/light` parameters as `GOTMLight` keywords `(; A, g1, g2)`;
- `initial_values`: every state variable's initialization, the Fortran default where the file
  gives none.

Only what this port implements is accepted: `couple_co2` must be true and `use_aice` false, and
an unknown parameter is an error.
"""
function fabm_configuration(path, FT::DataType = Float64)
    instances = YAML.load_file(path)["instances"]

    ecosmo = only_instance(instances, ECOSMO_MODEL, path)
    carbonate = only_instance(instances, CARBONATE_MODEL, path)
    light = only_instance(instances, LIGHT_MODEL, path)

    ecosmo_parameters = lowercase_keys(get(ecosmo, "parameters", Dict()))
    yaml_bool(get(ecosmo_parameters, "couple_co2", false)) ||
        throw(ArgumentError("$path: only `couple_co2: true` is implemented"))
    for instance_parameters in (ecosmo_parameters, lowercase_keys(get(carbonate, "parameters", Dict())))
        yaml_bool(get(instance_parameters, "use_aice", false)) &&
            throw(ArgumentError("$path: `use_aice: true` is not implemented"))
    end

    names = Dict(lowercase(String(name)) => name for name in ECOSMO_PARAMETER_NAMES)
    kwargs = Pair{Symbol, Float64}[]
    for (key, value) in ecosmo_parameters
        key in ("couple_co2", "use_aice") && continue
        haskey(names, key) || throw(ArgumentError("$path: unknown niva_ecosmo parameter `$key`"))
        push!(kwargs, names[key] => yaml_float(value))
    end
    parameters = EcosmoParameters(FT; kwargs...)

    carbonate_parameters = lowercase_keys(get(carbonate, "parameters", Dict()))
    pH_scale = Int(get(carbonate_parameters, "phscale", 1))

    light_parameters = lowercase_keys(get(light, "parameters", Dict()))
    light_kwargs = (A = yaml_float(get(light_parameters, "a", 0.58)),
                    g1 = yaml_float(get(light_parameters, "g1", 0.35)),
                    g2 = yaml_float(get(light_parameters, "g2", 23.0)))

    initialization = merge(lowercase_keys(get(ecosmo, "initialization", Dict())),
                           lowercase_keys(get(carbonate, "initialization", Dict())))
    initial_values = map(keys(ECOSMO_DEFAULT_INITIAL_VALUES), values(ECOSMO_DEFAULT_INITIAL_VALUES)) do name, default
        yaml_float(get(initialization, String(name), default))
    end

    return (; parameters, pH_scale, light = light_kwargs,
              initial_values = NamedTuple{keys(ECOSMO_DEFAULT_INITIAL_VALUES)}(initial_values))
end

function only_instance(instances, model, path)
    matches = [instance for instance in values(instances) if get(instance, "model", nothing) == model]
    length(matches) == 1 ||
        throw(ArgumentError("$path: expected one instance of `$model`, found $(length(matches))"))
    return matches[1]
end

lowercase_keys(dict) = Dict(lowercase(String(key)) => value for (key, value) in dict)

# YAML 1.1 reads `1e-7` (no decimal point) as a string; FABM reads it as a number.
yaml_float(value::Real) = Float64(value)
yaml_float(value::AbstractString) = parse(Float64, value)

yaml_bool(value::Bool) = value
yaml_bool(value::AbstractString) = lowercase(value) in ("true", ".true.", "yes", "on")
