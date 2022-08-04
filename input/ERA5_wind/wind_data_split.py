import xarray as xr
import matplotlib.pyplot as plt
ds=xr.open_dataset('ERA5_IO_01_91_21.nc')


d_split=[]
for i in range(1991,2021,5):
    dt=ds.sel(TIME=slice(repr(i)+'-01-01',repr(i+5)+'-01-01'))
    d_split.append(dt)
    
d_split[0].to_netcdf('era5_IO_01_91_95.nc')
d_split[1].to_netcdf('era5_IO_01_96_00.nc')
d_split[2].to_netcdf('era5_IO_01_01_05.nc')
d_split[3].to_netcdf('era5_IO_01_06_10.nc')
d_split[4].to_netcdf('era5_IO_01_11_15.nc')
d_split[5].to_netcdf('era5_IO_01_16_20.nc')
