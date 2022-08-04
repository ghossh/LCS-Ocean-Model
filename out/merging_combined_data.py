import xarray as xr
import glob
path="/home/NCAOR/supriyog/model_lcs/test_lcs_code0.1/out/out_files/**/combined_out.nc"
#path="/home/NCAOR/supriyog/model_lcs/test_lcs_code0.1/out/old_run1/**/combined_out.nc"
df=[]

# append datasets into the list
for i in glob.glob(path, recursive=True):
    print(i)
    temp_df = xr.open_mfdataset(i)
    df.append(temp_df)
    #print(i,"----------")


df_final= df[0]

for i in range(1,len(df)):
    #print(df[i])
    print(i,"-------------")
    df_final=df_final.merge(df[i])


df_final.to_netcdf("LCSCR_final_era5.nc")

