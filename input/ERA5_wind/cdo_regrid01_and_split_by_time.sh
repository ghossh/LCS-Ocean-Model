#cdo remapbil,./../era5_91_92.nc /home/NCAOR/supriyog/raw_data/ERA5_wind/ERA5_1991_21_daily.nc ERA5_IO_01_91_21.nc
cdo seldate,1991-01-01,1996-01-01 ERA5_IO_01_91_21.nc era5_IO_01_91_95.nc
cdo seldate,1996-01-01,2001-01-01 ERA5_IO_01_91_21.nc era5_IO_01_96_00.nc
cdo seldate,2001-01-01,2006-01-01 ERA5_IO_01_91_21.nc era5_IO_01_01_05.nc
cdo seldate,2006-01-01,2011-01-01 ERA5_IO_01_91_21.nc era5_IO_01_06_10.nc
cdo seldate,2011-01-01,2016-01-01 ERA5_IO_01_91_21.nc era5_IO_01_11_15.nc
cdo seldate,2016-01-01,2021-01-01 ERA5_IO_01_91_21.nc era5_IO_01_16_20.nc
