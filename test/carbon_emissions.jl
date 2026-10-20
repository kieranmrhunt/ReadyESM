const R=ReadyESM
const O=R.Oceananigans

function emissions_config(;rate=10.,mode=:emissions_driven,device=:cpu)
    R.ExperimentConfig(;device,ocean_model=:oceananigans,ocean_dynamics=true,
        forcing=R.ForcingConfig(radiation=:rrtmgp_clear_sky),
        carbon=R.CarbonConfig(;mode,emissions_pgc_per_year=rate))
end
function emissions_controller(config=emissions_config())
    ppm=config.carbon.atmospheric_co2_ppm
    R.CarbonController(config.carbon.mode,ppm*R._CARBON_PGC_PER_PPM,ppm,
        38000.,0.,0.,0.,0.,NaN,-1,900.,R._carbon_provenance(config))
end

@testset "Emissions mode configuration and offsets" begin
    for device in (:cpu,:gpu),rate in (0.,10.,1000.)
        cfg=emissions_config(;device,rate)
        @test R.validate(cfg;check_input_files=false)===cfg
        @test R._carbon_ocean_kwargs(cfg)==(tracers=(:T,:S,:DIC,:TA),)
        p=R._carbon_provenance(cfg)
        @test p.emissions_pgc_per_year==rate
        @test p.seconds_per_year==365.25*86400
        @test p.dic_storage_offset_umol_kg==2000.
        @test p.ta_storage_offset_umol_kg==2300.
    end
    for mode in (:off,:tracer_only,:prescribed)
        @test_throws ArgumentError R.validate(emissions_config(;mode);check_input_files=false)
    end
    for rate in (-1.,NaN,Inf)
        @test_throws ArgumentError R.validate(emissions_config(;rate);check_input_files=false)
    end
    @test R._carbon_storage_offset(:emissions_driven,Val(:DIC),Float32)==2000f0
    @test R._carbon_storage_offset(:emissions_driven,Val(:TA),Float32)==2300f0
    @test !hasproperty(R._carbon_provenance(emissions_config(mode=:prescribed,rate=0.)),:emissions_pgc_per_year)
end

@testset "Emissions ledger and signed exchange" begin
    for rate in (0.,10.,1000.),flux in (-1e-6,0.,1e-6)
        c=emissions_controller(emissions_config(;rate))
        atmosphere=c.atmospheric_carbon_pgc
        total=atmosphere+c.ocean_carbon_pgc
        R._advance_carbon_ledger!(c,0,0.,flux)
        for step in 1:96
            old=c.atmospheric_carbon_pgc
            R._advance_carbon_ledger!(c,step,900.0*step,flux)
            emitted=rate*(900.0*step/R._CARBON_EMISSIONS_SECONDS_PER_YEAR)
            @test isapprox(c.cumulative_emissions_pgc,emitted;rtol=32eps(Float64),atol=0)
            @test isapprox(c.atmospheric_carbon_pgc,atmosphere+emitted-flux*900*step;rtol=64eps(Float64))
            @test iszero(c.cumulative_implied_emissions_pgc)
            @test c.atmospheric_co2_ppm==c.atmospheric_carbon_pgc/R._CARBON_PGC_PER_PPM
            @test isapprox(c.atmospheric_carbon_pgc+c.ocean_carbon_pgc-c.cumulative_emissions_pgc,total;rtol=128eps(Float64))
            state=O.prognostic_state(c)
            R._advance_carbon_ledger!(c,step,900.0*step,99.)
            @test isequal(O.prognostic_state(c),state)
        end
    end
end

@testset "Rejected emissions ledger updates leave state unchanged" begin
    for (iteration,time,flux) in ((2,1800.,0.),(1,901.,0.),(1,900.,NaN),(0,1.,0.))
        c=emissions_controller()
        R._advance_carbon_ledger!(c,0,0.,0.)
        state=O.prognostic_state(c)
        @test_throws Exception R._advance_carbon_ledger!(c,iteration,time,flux)
        @test isequal(O.prognostic_state(c),state)
    end
    for flux in (1e9,-1e9)
        c=emissions_controller()
        R._advance_carbon_ledger!(c,0,0.,flux)
        state=O.prognostic_state(c)
        @test_throws ArgumentError R._advance_carbon_ledger!(c,1,900.,0.)
        @test isequal(O.prognostic_state(c),state)
    end
end

@testset "Emissions controller restart" begin
    c=emissions_controller()
    R._advance_carbon_ledger!(c,0,0.,1e-6)
    R._advance_carbon_ledger!(c,1,900.,-1e-6)
    state=O.prognostic_state(c)
    restored=emissions_controller()
    O.restore_prognostic_state!(restored,state)
    @test isequal(O.prognostic_state(restored),state)
    for step in 2:8
        flux=step*1e-7
        R._advance_carbon_ledger!(c,step,900.0*step,flux)
        R._advance_carbon_ledger!(restored,step,900.0*step,flux)
        @test isequal(O.prognostic_state(c),O.prognostic_state(restored))
    end
    @test_throws ErrorException O.restore_prognostic_state!(emissions_controller(emissions_config(rate=11.)),state)
    @test_throws ErrorException O.restore_prognostic_state!(emissions_controller(),merge(state,(cumulative_implied_emissions_pgc=1.,)))
    @test_throws ErrorException O.restore_prognostic_state!(emissions_controller(),merge(state,(cumulative_emissions_pgc=-1.,)))
end
