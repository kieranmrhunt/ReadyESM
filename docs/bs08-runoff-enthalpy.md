# Runoff enthalpy (BS08, optional)

## Known problem

Keep `terrarium_runoff_enthalpy: false` in climate runs. In two matched
31-day tests, enabling the option raised subsurface soil temperatures to
192 °C. Both controls stayed below 43 °C. The river heat budget balanced,
but heat transport inside the soil was incomplete. A repair is under test.

## What it does

By default, land runoff reaches the ocean at the sea surface temperature of the receiving cell. The ocean
receives the water, but the heat of the water comes from the ocean itself.

With `terrarium_runoff_enthalpy: true`, the land model keeps the heat content of its surface water. The heat moves
with the water from the soil surface, through the runoff store, to the ocean river mouths. The ocean then receives
the heat that the water actually carries. The reference state is liquid water at 0 °C.

The option is off by default. With the option off, the model and its output are unchanged.

## How to use it

For isolated development tests, enable the option in the experiment configuration:

```yaml
terrarium_runoff_enthalpy: true
```

The option needs the Terrarium land model and the mutable z-star ocean grid. The model stops at construction if
these are not present.

## Output

With the option on, the NetCDF diagnostics contain these variables. Their names start with `bs08_`.

| Variable | Units | Meaning |
| --- | --- | --- |
| `bs08_source_enthalpy_w` | W | Runoff enthalpy rate that leaves the land |
| `bs08_receiver_enthalpy_w` | W | Runoff enthalpy rate at the ocean receivers |
| `bs08_source_water_kgs`, `bs08_receiver_water_kgs` | kg s⁻¹ | Routed water at the land and at the ocean |
| `bs08_routed_water_mean_temperature_c` | °C | Mean temperature of the routed water. It is missing when no water flows. |
| `bs08_zero_flow` | 1 | 1 when no routed water flows |
| `bs08_minus_native_heat_w` | W | Receiver enthalpy minus the default heat term |
| `bs08_cumulative_released_enthalpy_j` | J | Enthalpy released from the land since launch |
| `bs08_cumulative_received_enthalpy_j` | J | Enthalpy applied to the ocean since launch |
| `bs08_surface_excess_enthalpy_j`, `bs08_pending_runoff_enthalpy_j` | J | Enthalpy in the land surface stores |
| `bs08_cumulative_clamp_enthalpy_j` | J | Enthalpy given to water that the bottom saturation limit creates |
| `bs08_partition_inconsistent_columns` | 1 | Land columns where the runoff law is not the default law. It must be 0. The launch row comes before the first land step and shows 0. |
| `bs08_receiver_enthalpy_flux` | W m⁻² | Final receiver enthalpy flux on the ocean grid |

The time series use the same sample times as the coupled budget series.

The ocean applies the routed enthalpy one coupled step after the land releases it. Thus the released enthalpy
minus the received enthalpy is the enthalpy in transit. The file attribute `bs08_in_transit_enthalpy_j` gives this
value at the end of the run.

This attribute describes the run's current state. It can differ between
forcing experiments because the runoff heat changes.

## Restart

The runoff enthalpy stores and the diagnostic history are in the restart file. A restart continues the history
without a gap. An older restart file without the history can start a run only at launch. A run that continues an
older file stops at its first diagnostic sample.

## Qualification

CPU comparisons check unchanged off-mode behaviour, routed water and heat,
output histories and exact restarts. The full-grid A100 test checks the
active energy ledger and GPU storage. The joint test also enables prescribed
methane and nitrous oxide. See the [validation record](validation-runoff-trace-gases.md).

An earlier four-block coupled timing test measured a mean on/off ratio of
1.0044. A later short coupled pair measured 1.076. The land-only kernel test
exceeded its 1.25 limit. These tests used different harnesses and source
revisions. Allow for extra runtime when enabling the option.

The monthly comparison found the soil-temperature problem described above.
The option remains off by default.
