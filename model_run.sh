source /home/NCAOR/supriyog/miniconda3/bin/activate
conda activate clim2

gitpush
rm -r nohup.out
cp -dir ./header_files ./header_files_1
rm -r ./restart/*
rm -r ./out/*.nc
rm -r ./header_files_before_restart
 cp -dir ./header_files ./header_files_before_restart
rm -rf ./out/out_files 
mkdir ./out/out_files
rm -rf ./input/*restart*
rm -rf ./input/restart_files 
mkdir ./input/restart_files

# gitpush
# ## 1st run ##
./submit
python ./out/file_merge.py
cp ./restart/* ./input
mkdir ./input/restart_files/run1
cp ./restart/* ./input/restart_files/run1
mkdir ./out/out_files/run1
mv ./out/*nc ./out/out_files/run1
#gitpush

# ## 2nd run ##
rm -r ./header_files
cp -dir ./header_files_before_restart ./header_files
python changes_for_longrun.py era5_IO_01_96_00.nc 1 1826
./submit
rm ./input/*restart*
python ./out/file_merge.py
cp ./restart/* ./input
mkdir ./input/restart_files/run2
cp ./restart/* ./input/restart_files/run2
mkdir ./out/out_files/run2
mv ./out/*nc ./out/out_files/run2
#gitpush

# 3rd run ##
rm -r ./header_files
cp -dir ./header_files_before_restart ./header_files
python changes_for_longrun.py era5_IO_01_01_05.nc 1 1826
./submit
rm ./input/*restart*
python ./out/file_merge.py
cp ./restart/* ./input
mkdir ./input/restart_files/run3
cp ./restart/* ./input/restart_files/run3
mkdir ./out/out_files/run3
mv ./out/*nc ./out/out_files/run3
gitpush

# 4th run ##
rm -r ./header_files
cp -dir ./header_files_before_restart ./header_files
python changes_for_longrun.py era5_IO_01_06_10.nc 1 1826
./submit
rm ./input/*restart*
python ./out/file_merge.py
cp ./restart/* ./input
mkdir ./input/restart_files/run4
cp ./restart/* ./input/restart_files/run4
mkdir ./out/out_files/run4
mv ./out/*nc ./out/out_files/run4
#gitpush

## 5th run ##
rm -r ./header_files
cp -dir ./header_files_before_restart ./header_files
python changes_for_longrun.py era5_IO_01_11_15.nc 1 1826
./submit
rm ./input/*restart*
python ./out/file_merge.py
cp ./restart/* ./input
mkdir ./input/restart_files/run5
cp ./restart/* ./input/restart_files/run5
mkdir ./out/out_files/run5
mv ./out/*nc ./out/out_files/run5
#gitpush

## 6th run ##
rm -r ./header_files
cp -dir ./header_files_before_restart ./header_files
python changes_for_longrun.py era5_IO_01_16_20.nc 1 1827
./submit
rm ./input/*restart*
python ./out/file_merge.py
cp ./restart/* ./input
mkdir ./input/restart_files/run6
cp ./restart/* ./input/restart_files/run6
mkdir ./out/out_files/run6
mv ./out/*nc ./out/out_files/run6
#gitpush




# reset
rm -r ./header_files
rm ./input/*restart*
cp -dir ./header_files_before_restart ./header_files
cd ./out
python merging_combined_data.py
#gitpush
