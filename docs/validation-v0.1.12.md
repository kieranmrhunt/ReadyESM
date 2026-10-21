# Optional moist-mass coupling

The new `atmosphere_moist_mass_coupling` switch lets evaporation and
precipitation change atmospheric mass. It supports the Betts–Miller,
net-column-precipitation and RRTMGP host. The switch defaults to `false`.

In one paired T31/L27 January run, dry-air change over 31 days was
−6.51 kg/m² with the switch off and −0.24 kg/m² with it on. Column-water
RMSE moved from 5.989 to 5.959 kg/m². T500 RMSE moved from 4.770 to
4.777 K. Precipitation RMSE moved from 2.492 to 2.479 mm/day. Net TOA
changed by −0.38 W/m².

The v0.1.11 port passed the full CPU suite with 54 testsets, against
49 in the control. The off state and output match v0.1.11 at every
checked stage. An on-state fresh-process restart is exact. Mismatched
old and new checkpoints are rejected before state changes. On an A100,
15/15 checks pass, and the eight-step dry-air change is −5.06 kg/m²
off versus +0.096 kg/m² on. The column source matches evaporation minus
precipitation to 5.8e−7 relative in the short coupled check.

The on state retains about −0.2 kg/m² dry-air change during launch.
Water-source enthalpy is recorded as a diagnostic and is not applied
to the energy equations. The option stays off in production while
longer coupled skill and energy accounting are tested.
