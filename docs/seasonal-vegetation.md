# Seasonal ERA5 vegetation cover and leaf area

The v0.1.19 production profile selects `terrarium_vegetation_cycle: cover_lai`
for prescribed seasonal vegetation. The generic experiment default `off`
retains the configured vegetation inputs. The cycle requires Terrarium, ERA5
instantaneous land initialization, and
`terrarium_evapotranspiration: era5_prescribed_vegetation`.

Use the versioned 1993 monthly land-cycle file and its recorded SHA-256:

```yaml
terrarium_initial_conditions: era5_instantaneous
terrarium_evapotranspiration: era5_prescribed_vegetation
terrarium_vegetation_cycle: cover_lai
terrarium_vegetation_cycle_path: data/era5_seasonal_land_cycle_1993/era5_monthly_land_cycle_1993.nc
terrarium_vegetation_cycle_sha256: b80d7e01a4390340794030f1f924e5583d813fb163a1f9d875ae749091c13f60
```

The monthly fields are remapped to the model's native land columns and interpolated between Gregorian month centres. ReadyESM updates vegetation cover and leaf area index at each Terrarium timestep. It preserves root fractions and land prognostic stores. The source and native mapping identities are checked on restart.

The production profile also enables the separate [seasonal albedo](seasonal-albedo.md)
and dissipative-heating settings. Its combined 90-day run cannot attribute a
temperature change to vegetation alone; separate paired arms are needed to
assess individual effects.
