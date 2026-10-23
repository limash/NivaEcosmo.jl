using NivaEcosmo
using NivaEcosmo: interior_rates, bottom_fluxes, oxygen_surface_flux, carbonate_diagnostics,
    carbon_dioxide_surface_flux, day_length, gotm_light, attenuation_coefficient, ECOSMO_PARAMETER_NAMES,
    PER_DAY, MGC_PER_MMOL_N, carbonate_system
using Oceananigans
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
        rates = interior_rates(p, c, 10.0, 0.6, 50.0, 2.5).rates
        gross = sum(name -> abs(rates[name]), P_TRACERS)
        @test abs(sum(name -> rates[name], P_TRACERS)) < 1e-12 * gross
    end

    @testset "carbonate iteration" begin
        cs = carbonate_system(10.0, 30.0, 1.0, 2.0e-3, 2.2e-3, 1)
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
        D = day_length(latitude, yday)
        ext = attenuation_coefficient(p, state.flachl, state.diachl)
        light = gotm_light(swr, configuration.light.A, configuration.light.g1, configuration.light.g2,
                           0.0 + 10.0 * 0.5, 0.0 + ext * 10.0 * 0.5)
        previous = (pH = 0.0, pCO2 = 0.0, CarbA = 0.0, BiCarb = 0.0, Carb = 0.0)
        carbonate = carbonate_diagnostics(state.dic, state.alk, T, S, ρ, 5.0, previous, configuration.pH_scale)
        interior = interior_rates(p, state, T, D, light.par, carbonate.Om_cal)
        surface = (oxy = oxygen_surface_flux(state.oxy, T, S, wind),
                   dic = carbon_dioxide_surface_flux(state.dic, state.alk, T, S, ρ, 5.0, wind, 420.0,
                                                     carbonate.pCO2, configuration.pH_scale))
        bottom = bottom_fluxes(p, state, state, T, carbonate.Om_cal, 0.0)

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
        # Sinking through an immersed seabed into the sediment, with burial switched off: the
        # phosphorus in the water plus the sediment is constant to round-off.
        p = EcosmoParameters(; burialRt = 0, burialRt2 = 0, burialRtO = 0, burialRt2O = 0, burialRtC = 0,
                               burialRt2C = 0, cz_sinkD = 0.5)
        underlying = LatitudeLongitudeGrid(size = (4, 3, 6), longitude = (10, 10.4), latitude = (59, 59.3),
                                           z = (-30, 0), topology = (Bounded, Bounded, Bounded))
        grid = ImmersedBoundaryGrid(underlying, GridFittedBottom((λ, φ) -> -12 - 40 * (λ - 10)))

        surface_shortwave(λ, φ, t) = 200.0
        bgc = ECOSMO(grid; parameters = p, surface_shortwave, start_day_of_year = 150,
                     sediment = ecosmo_sediment(grid; parameters = p, bottom_stress = QuadraticBottomStress()))
        boundary_conditions = map(bc -> FieldBoundaryConditions(top = bc), ecosmo_surface_fluxes(; wind_speed = 5.0))
        model = HydrostaticFreeSurfaceModel(grid; tracers = (:T, :S), biogeochemistry = bgc, boundary_conditions)

        init = fabm_configuration(BOX_YAML).initial_values
        set!(model; T = 10, S = 30, (name => init[name] for name in ECOSMO_TRACERS)...)
        set!(model; det = 50.0)
        for name in (:sed1, :sed2, :sed3, :sed4)
            set!(bgc.sediment.fields[name], init[name])
        end

        function total_phosphorus(model)
            Nx, Ny, Nz = size(model.grid)
            water = sum(i -> sum(j -> sum(k -> immersed_cell(i, j, k, model.grid) ? 0.0 :
                Vᶜᶜᶜ(i, j, k, model.grid) * sum(name -> model.tracers[name][i, j, k], P_TRACERS), 1:Nz), 1:Ny), 1:Nx)
            sediment = model.biogeochemistry.sediment.fields
            bottom = sum(i -> sum(j -> Azᶜᶜᶜ(i, j, 1, model.grid) * (sediment.sed1[i, j, 1] + sediment.sed3[i, j, 1]), 1:Ny), 1:Nx)
            return water + bottom
        end

        P₀ = total_phosphorus(model)
        sed1₀ = sum(interior(bgc.sediment.fields.sed1))
        for _ in 1:20
            time_step!(model, 600)
        end

        @test all(name -> all(isfinite, interior(model.tracers[name])), ECOSMO_TRACERS)
        @test sum(interior(bgc.sediment.fields.sed1)) > sed1₀  # detritus has settled
        @test total_phosphorus(model) ≈ P₀ rtol = 1e-12
    end
end
