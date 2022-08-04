import xarray as xr
import glob
path="./out/*out.nc"
df=[]

# append datasets into the list
for i in glob.glob(path, recursive=True):
    temp_df = xr.open_mfdataset(i)
    df.append(temp_df)
    #print(i,"----------")
  

df_final= 0

for i in range(len(df)):
   # print(df[i])
   # print(i,"-------------")
    df_final=df_final+df[i]
    

df_final.to_netcdf("out/combined_out.nc")
