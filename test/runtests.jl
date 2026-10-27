using NivaEcosmo
using NivaEcosmo: interior_rates, bottom_fluxes, oxygen_surface_flux, carbonate_diagnostics,
    carbon_dioxide_surface_flux, day_length, gotm_light, attenuation_coefficient, ECOSMO_PARAMETER_NAMES,
    PER_DAY, MGC_PER_MMOL_N, carbonate_system, ecosmo_constants
using Oceananigans
using Oceananigans.BoundaryConditions: getbc
using Oceananigans.TimeSteppers: update_state!
using Oceananigans.Models.HydrostaticFreeSurfaceModels: compute_tracer_tendencies!
using Oceananigans.Operators: Vᶜᶜᶜ, Azᶜᶜᶜ
using Oceananigans.ImmersedBoundaries: immersed_cell
using Test

const BOX_YAML = joinpath(@__DIR__, "..", "data", "fabm.yaml")

# The box model's density, so a model cell and the pointwise functions see the same one.
box_density(T, S, P) = 1000.0 + 0.78 * S - 0.0065 * (T - 4.0)^2

# The phosphorus-bearing pools, all in mgC in Redfield units.
const P_TRACERS = (:pho, :fla, :dia, :microzoo, :mesozoo, :det, :dom)

# Sinking through an immersed seabed into the sediment, with resuspension and gas exchange and with
# burial switched off, on a `FT` grid. The shortwave function returns `Float64` whatever `FT` is.
function phosphorus_model(FT)
    p = EcosmoParameters(FT; burialRt = 0, burialRt2 = 0, burialRtO = 0, burialRt2O = 0, burialRtC = 0,
                             burialRt2C = 0, cz_sinkD = 0.5)
    underlying = LatitudeLongitudeGrid(FT; size = (4, 3, 6), longitude = (10, 10.4), latitude = (59, 59.3),
                                       z = (-30, 0), topology = (Bounded, Bounded, Bounded))
    grid = ImmersedBoundaryGrid(underlying, GridFittedBottom((λ, φ) -> -12 - 40 * (λ - 10)))

    surface_shortwave(λ, φ, t) = 200.0
    bgc = ECOSMO(grid; parameters = p, surface_shortwave, start_day_of_year = 150,
                 sediment = ecosmo_sediment(grid; parameters = p, bottom_stress = QuadraticBottomStress()))
    boundary_conditions = map(bc -> FieldBoundaryConditions(top = bc), ecosmo_surface_fluxes(FT; wind_speed = 5.0))
    model = HydrostaticFreeSurfaceModel(grid; tracers = (:T, :S), biogeochemistry = bgc, boundary_conditions)

    init = fabm_configuration(BOX_YAML).initial_values
    set!(model; T = 10, S = 30, (name => init[name] for name in ECOSMO_TRACERS)...)
    set!(model; det = 50.0)
    for name in (:sed1, :sed2, :sed3, :sed4)
        set!(bgc.sediment.fields[name], init[name])
    end

    return model
end

# Water-column plus sediment phosphorus, summed in `Float64` whatever the model's type.
function total_phosphorus(model)
    Nx, Ny, Nz = size(model.grid)
    water = sum(i -> sum(j -> sum(k -> immersed_cell(i, j, k, model.grid) ? 0.0 :
        Vᶜᶜᶜ(i, j, k, model.grid) * sum(name -> Float64(model.tracers[name][i, j, k]), P_TRACERS),
        1:Nz), 1:Ny), 1:Nx)
    sediment = model.biogeochemistry.sediment.fields
    bottom = sum(i -> sum(j -> Azᶜᶜᶜ(i, j, 1, model.grid) *
        (Float64(sediment.sed1[i, j, 1]) + Float64(sediment.sed3[i, j, 1])), 1:Ny), 1:Nx)
    return water + bottom
end

@testset "NivaEcosmo" begin
    @testset "parameters" begin
        p = EcosmoParameters()
        @test p.muPs == 1.10 * PER_DAY
        @test p.rNO3Ps == 0.50 * MGC_PER_MMOL_N
        @test p.Exphy == 0.04
        @test EcosmoParameters(; muPs = 2.0).muPs == 2.0 * PER_DAY
        @test_throws ArgumentError EcosmoParameters(; muPS = 2.0)
    end

    @testset "fabm.yaml" begin
        configuration = fabm_configuration(BOX_YAML)
        @test configuration.parameters.muPs == 1.47 * PER_DAY
        @test configuration.parameters.burialRt2 == 1e-7 * PER_DAY  # written `1e-7`, which YAML reads as a string
        @test configuration.pH_scale == 1
        @test configuration.light == (A = 0.58, g1 = 0.35, g2 = 4.0)
        @test configuration.initial_values.no3 == 955.0
        @test configuration.initial_values.dic == 1450.0

        # ECOSMO takes its float type from the grid, so a Float32 grid needs Float32 parameters.
        grid = RectilinearGrid(Float32; size = (1, 1, 2), extent = (1, 1, 10), topology = (Periodic, Periodic, Bounded))
        parameters = fabm_configuration(BOX_YAML, Float32).parameters
        @test ECOSMO(grid; parameters, latitude = 60).underlying_biogeochemistry.parameters === parameters
    end

    @testset "interior phosphorus is conserved" begin
        p = fabm_configuration(BOX_YAML).parameters
        c = (no3 = 900.0, nh4 = 80.0, pho = 700.0, sil = 1900.0, oxy = 300.0, fla = 40.0, dia = 60.0,
             flachl = 0.8, diachl = 1.5, microzoo = 10.0, mesozoo = 8.0, det = 30.0, opa = 20.0, cal = 1.0,
             dom = 50.0)
        rates = interior_rates(p, ecosmo_constants(), c, 10.0, 0.6, 50.0, 2.5).rates
        gross = sum(name -> abs(rates[name]), P_TRACERS)
        @test abs(sum(name -> rates[name], P_TRACERS)) < 1e-12 * gross
    end

    @testset "carbonate iteration" begin
        cs = carbonate_system(ecosmo_constants(), 10.0, 30.0, 1.0, 2.0e-3, 2.2e-3, 1)
        @test cs.success
        @test 7.5 < cs.pH < 8.5
    end

    @testset "Oceananigans tendencies equal the pointwise rates" begin
        # One 10 m cell is the box model: its tendencies must be the interior rates plus the
        # surface and bottom fluxes over the cell thickness.
        configuration = fabm_configuration(BOX_YAML)
        p = configuration.parameters
        init = configuration.initial_values
        T, S, swr, wind, yday, latitude = 8.0, 25.0, 150.0, 6.0, 100.25, 59.7

        grid = RectilinearGrid(size = 1, z = (-10, 0), topology = (Flat, Flat, Bounded))
        bgc = ECOSMO(grid; parameters = p, light_attenuation = GOTMLight(grid, swr; configuration.light...),
                     sediment = ecosmo_sediment(grid; parameters = p), latitude, start_day_of_year = yday,
                     pH_scale = configuration.pH_scale, density = box_density)
        surface_fluxes = ecosmo_surface_fluxes(; wind_speed = wind, pH_scale = configuration.pH_scale, density = box_density)
        boundary_conditions = map(bc -> FieldBoundaryConditions(top = bc), surface_fluxes)
        model = HydrostaticFreeSurfaceModel(grid; tracers = (:T, :S), biogeochemistry = bgc, boundary_conditions)

        set!(model; T, S, (name => init[name] for name in ECOSMO_TRACERS)...)
        sediment = bgc.sediment
        for name in (:sed1, :sed2, :sed3, :sed4)
            set!(sediment.fields[name], init[name])
        end

        update_state!(model)
        compute_tracer_tendencies!(model)

        # the box model's composition of the pointwise functions
        state = init
        ρ = box_density(T, S, 0.5)
        constants = ecosmo_constants()
        D = day_length(constants, latitude, yday)
        ext = attenuation_coefficient(p, state.flachl, state.diachl)
        light = gotm_light(swr, configuration.light.A, configuration.light.g1, configuration.light.g2,
                           0.0 + 10.0 * 0.5, 0.0 + ext * 10.0 * 0.5)
        previous = (pH = 0.0, pCO2 = 0.0, CarbA = 0.0, BiCarb = 0.0, Carb = 0.0)
        carbonate = carbonate_diagnostics(constants, state.dic, state.alk, T, S, ρ, 5.0, previous,
                                          configuration.pH_scale)
        interior = interior_rates(p, constants, state, T, D, light.par, carbonate.Om_cal)
        surface = (oxy = oxygen_surface_flux(constants, state.oxy, T, S, wind),
                   dic = carbon_dioxide_surface_flux(constants, state.dic, state.alk, T, S, ρ, 5.0, wind, 420.0,
                                                     carbonate.pCO2, configuration.pH_scale))
        bottom = bottom_fluxes(p, constants, state, state, T, carbonate.Om_cal, 0.0)

        @test bgc.underlying_biogeochemistry.day_length[1, 1, 1] == D
        @test Oceananigans.fields(model).PAR[1, 1, 1] == light.par
        @test bgc.underlying_biogeochemistry.carbonate.Om_cal[1, 1, 1] == carbonate.Om_cal

        for name in ECOSMO_TRACERS
            expected = interior.rates[name] + get(surface, name, 0.0) / 10.0 + get(bottom.fluxes, name, 0.0) / 10.0
            @test model.timestepper.Gⁿ[name][1, 1, 1] ≈ expected rtol = 1e-12
        end
        for name in (:sed1, :sed2, :sed3, :sed4)
            @test sediment.timestepper.Gⁿ[name][1, 1] == bottom.sediment[name]
        end
    end

    @testset "3D run conserves phosphorus" begin
        # The phosphorus in the water plus the sediment is constant to round-off.
        model = phosphorus_model(Float64)
        bgc = model.biogeochemistry

        P₀ = total_phosphorus(model)
        sed1₀ = sum(interior(bgc.sediment.fields.sed1))
        for _ in 1:20
            time_step!(model, 600)
        end

        @test all(name -> all(isfinite, interior(model.tracers[name])), ECOSMO_TRACERS)
        @test sum(interior(bgc.sediment.fields.sed1)) > sed1₀  # detritus has settled
        @test total_phosphorus(model) ≈ P₀ rtol = 1e-12
    end

    @testset "Float32 computes in Float32" begin
        # A stray `Float64` number would promote these results.
        F = Float32
        p = fabm_configuration(BOX_YAML, F).parameters
        constants = ecosmo_constants(F)
        init = fabm_configuration(BOX_YAML).initial_values
        reactive = NivaEcosmo.ECOSMO_REACTIVE_TRACERS
        c = map(F, NamedTuple{reactive}(map(name -> init[name], reactive)))
        sediment = map(F, (sed1 = init.sed1, sed2 = init.sed2, sed3 = init.sed3, sed4 = init.sed4))
        previous = map(F, (pH = 8.0, pCO2 = 400.0, CarbA = 0.0, BiCarb = 0.0, Carb = 100.0))

        @test all(x -> x isa F, interior_rates(p, constants, c, F(8), F(0.6), F(50), F(2.5)).rates)
        @test all(x -> x isa F, interior_rates(p, constants, c, F(8), F(0.6), F(0), F(2.5)).rates)
        @test all(x -> x isa F, bottom_fluxes(p, constants, c, sediment, F(8), F(2.5), F(0.1)).fluxes)
        for pH_scale in (-1, 0, 1)
            @test all(x -> x isa F, carbonate_diagnostics(constants, F(1450), F(1600), F(8), F(25), F(1020), F(5),
                                                          previous, pH_scale))
        end
        @test oxygen_surface_flux(constants, F(300), F(8), F(25), F(6)) isa F
        @test carbon_dioxide_surface_flux(constants, F(1450), F(1600), F(8), F(25), F(1020), F(5), F(6), F(420),
                                          F(400), 1) isa F
        @test day_length(constants, F(59.7), 100.25) isa F
        @test gotm_light(F(150), F(0.58), F(0.35), F(4), F(5), F(0.1)).par isa F

        model = phosphorus_model(F)
        grid = model.grid
        bgc = model.biogeochemistry.underlying_biogeochemistry
        model_fields = Oceananigans.fields(model)
        @test bgc(2, 2, 6, grid, Val(:no3), model.clock, model_fields, model_fields) isa F
        @test getbc(model.tracers.oxy.boundary_conditions.top, 2, 2, grid, model.clock, model_fields) isa F
        @test getbc(model.tracers.dic.boundary_conditions.top, 2, 2, grid, model.clock, model_fields) isa F
    end

    @testset "a Float32 run follows the Float64 one" begin
        # Twenty 10-minute steps of the conservation test's model in both precisions. Measured: every
        # tracer and sediment pool within 7e-7 of its largest value, phosphorus conserved to 3e-8.
        model64 = phosphorus_model(Float64)
        model32 = phosphorus_model(Float32)
        P₀ = total_phosphorus(model32)
        for _ in 1:20
            time_step!(model64, 600)
            time_step!(model32, 600)
        end

        sediment64 = model64.biogeochemistry.sediment.fields
        sediment32 = model32.biogeochemistry.sediment.fields
        fields64 = [(model64.tracers[name] for name in ECOSMO_TRACERS)..., values(sediment64)...]
        fields32 = [(model32.tracers[name] for name in ECOSMO_TRACERS)..., values(sediment32)...]
        for (a, b) in zip(fields64, fields32)
            @test eltype(b) == Float32
            @test maximum(abs, interior(a) .- interior(b)) <= 1e-5 * maximum(abs, interior(a))
        end
        @test total_phosphorus(model32) ≈ P₀ rtol = 1e-6
    end
end
