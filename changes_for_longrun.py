import os,glob,sys

##### ERA5 data wind ####
ch1=   """windfile='./input/era5_IO_01_91_95.nc'""" 
ch1n=  """windfile='./input/"""+sys.argv[1]+"""'"""


for f in glob.glob("./header_files/input_files.h", recursive=True):	
	print(f)
	with open(f) as file:
		s = file.read()
		s = s.replace(ch1,ch1n)
	with open(f, 'w') as file:
		file.write(s)


ch3= """start_program=0""" 
ch3n= """start_program="""+sys.argv[2]

ch4="""ndrun=1826"""
ch4n="""ndrun="""+sys.argv[3]
# 730 1095
for f in glob.glob("./header_files/param.h", recursive=True):	
	print(f)
	with open(f) as file:
		s = file.read()
		s = s.replace(ch3,ch3n)
		s = s.replace(ch4,ch4n)

	with open(f, 'w') as file:
		file.write(s)


