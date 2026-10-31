"""Retrieve one annual cycle of ERA5 vegetation and snow-free albedo inputs.

Keep all four spectral/direct-diffuse albedos. A broadband conversion has
not been chosen. Forecast albedo is deliberately excluded because it
includes the model's snow/ice modifications.
"""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np
from netCDF4 import Dataset, num2date

HERE = Path(__file__).resolve().parent
OUT = Path(__file__).resolve().parent
DATASET = 'reanalysis-era5-single-levels-monthly-means'
VARIABLES = {
    'cvh': 'high_vegetation_cover',
    'cvl': 'low_vegetation_cover',
    'lai_hv': 'leaf_area_index_high_vegetation',
    'lai_lv': 'leaf_area_index_low_vegetation',
    'aluvp': 'uv_visible_albedo_for_direct_radiation',
    'aluvd': 'uv_visible_albedo_for_diffuse_radiation',
    'alnip': 'near_ir_albedo_for_direct_radiation',
    'alnid': 'near_ir_albedo_for_diffuse_radiation',
    'lsm': 'land_sea_mask',
}
REQUEST = dict(product_type=['monthly_averaged_reanalysis'],year=['1993'],
    month=[f'{m:02d}' for m in range(1,13)],time=['00:00'],
    area=[90,0,-90,359],grid=[1,1],variable=list(VARIABLES.values()),
    data_format='netcdf',download_format='unarchived')

def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()

def validate(path):
    with Dataset(path) as ds:
        tn=next(n for n in ('valid_time','time') if n in ds.variables)
        dates=num2date(ds[tn][:],ds[tn].units,calendar=getattr(ds[tn],'calendar','standard'))
        if [(t.year,t.month,t.day) for t in dates] != [(1993,m,1) for m in range(1,13)]:
            raise ValueError('expected all twelve 1993 monthly records in order')
        if not np.array_equal(ds['latitude'][:],np.arange(90,-90.1,-1)):
            raise ValueError('unexpected latitude grid')
        if not np.array_equal(ds['longitude'][:],np.arange(360)):
            raise ValueError('unexpected longitude grid or duplicate seam')
        fields={}
        metadata={}
        for name in VARIABLES:
            var=ds[name]
            if var.dimensions != (tn,'latitude','longitude') or var.shape != (12,181,360):
                raise ValueError(f'unexpected dimensions: {name}')
            values=np.ma.asarray(var[:],dtype=float).filled(np.nan)
            valid=np.isfinite(values)
            if not valid.any():
                raise ValueError(f'no finite values: {name}')
            units=str(getattr(var,'units','')).strip()
            allowed={'m**2 m**-2','m2 m-2','m^2 m^-2','m2/m2','1'} if name.startswith('lai_') else {'1','(0 - 1)','(0-1)','dimensionless','~'}
            if units not in allowed:
                raise ValueError(f'unexpected units for {name}: {units!r}')
            upper=10 if name.startswith('lai_') else 1
            if values[valid].min() < 0 or values[valid].max() > upper:
                raise ValueError(f'out-of-range values: {name}')
            fields[name]=values
            metadata[name]=dict(units=units,shape=list(values.shape),
                finite_fraction=float(valid.mean()),minimum=float(values[valid].min()),
                maximum=float(values[valid].max()))
        if not np.isfinite(fields['lsm']).all():
            raise ValueError('missing land-sea mask')
        land=fields['lsm']>0
        for name,values in fields.items():
            if not np.isfinite(values[land]).all():
                raise ValueError(f'missing source values over land: {name}')
        return metadata

def main():
    import cdsapi
    OUT.mkdir(parents=True,exist_ok=True)
    target=OUT/'era5_monthly_land_cycle_1993.nc'
    intent=OUT/'request_intent.json'
    record=OUT/'request.json'
    receipt=OUT/'receipt.json'
    if target.exists() or intent.exists() or record.exists() or receipt.exists():
        raise FileExistsError('existing request or data: inspect its remote ID before resuming')
    with intent.open('x') as stream:
        json.dump(dict(dataset=DATASET,request=REQUEST,script_sha256=digest(__file__),
            created_utc=datetime.now(timezone.utc).isoformat()),stream,indent=2)
    client=cdsapi.Client(quiet=True,wait_until_complete=False)
    remote=client.retrieve(DATASET,REQUEST)
    with record.open('x') as stream:
        json.dump(dict(dataset=DATASET,request=REQUEST,request_id=remote.request_id),stream,indent=2)
    print('SEASONAL_LAND_REQUEST',remote.request_id,flush=True)
    temporary=target.with_suffix('.nc.part')
    remote.download(str(temporary))
    metadata=validate(temporary)
    temporary.replace(target)
    with receipt.open('x') as stream:
        json.dump(dict(schema='era5_seasonal_land_source_v1',path=str(target),
            sha256=digest(target),script_sha256=digest(__file__),
            request_id=remote.request_id,variables=metadata,
            spectral_albedo_conversion='not_yet_selected'),stream,indent=2)
    print('SEASONAL_LAND_DOWNLOAD_PASS',target,digest(target),flush=True)

if __name__=='__main__':
    main()
