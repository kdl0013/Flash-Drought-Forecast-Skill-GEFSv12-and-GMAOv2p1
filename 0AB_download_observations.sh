#!/bin/bash

#Download SMERGE soil moisture, gridMET PET, EDDI evaporative demand, and elevation
#To download SubX data, use script download_SubX_4_models.sh (run in HPC)

#My own environment is spyder
cas 
module load cdo
module load nco

#MUST SET inputs: main_directory and number of processors for CPU parallelism

#You can make main_directory path any path you would like and this will get
#all other directories in order. 

main_directory='/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
processors=8

data_d=$main_directory/Data
data_s=$main_directory/Scripts
mkdir -p $data_s 

mask=$data_s/CONUS_mask
merra_dir=$data_d/MERRA2
############## MERRA 2 ##################################################
mkdir -p $data_d/MERRA2
python3 $data_s/00_download_MERRA2.py #make files to download with multiprocessing

for dir in $merra_dir/*;do
cd $dir 
find . -size 0 -delete;
done

bash $data_s/wget_MERRA.sh
#Sometimes, GES DISC doesn't have a good connection and makes empty files

bash $data_s/wget_MERRA_small.sh  #last 30 days download


ulimit -n 9000
#Merge files into 1
cdo -remapcon,$mask/CONUS_mask.grd -dayavg -mergetime $merra_dir/RZSM/*.nc4 $merra_dir/RZSM.nc4
cdo -remapcon,$mask/CONUS_mask.grd -dayavg -mergetime $merra_dir/radiation/*.nc4 $merra_dir/radiation.nc4
cdo -remapcon,$mask/CONUS_mask.grd -mergetime $merra_dir/temperature_min_mean_max/*.nc4 $merra_dir/temperature.nc4
cdo -remapcon,$mask/CONUS_mask.grd -dayavg -mergetime $merra_dir/wind_humidity/*.nc4 $merra_dir/wind_humidity.nc4
cdo -remapcon,$mask/CONUS_mask.grd -dayavg -mergetime $merra_dir/surface_RZSM/*.nc4 $merra_dir/surface_RZSM.nc4
############## Evaporative Demand Drought Index ##########################
{
echo "Working on EDDI pre-processing"

eddi=$data_d/EDDI

eddi_save_dir=${eddi}/convert_2_nc
mkdir -p $eddi_save_dir

cd ${eddi} #Move to directory so files save properly
bash $data_s/0c_download_EDDI_wget.sh
}
#Some downloads produce double on my computer, this moves the doubles so you
#could verify if the actual file is present

#for f in *.asc; 
#do  
#    if [[ ${#f} -gt 28 ]]; 
#then 
#    rm "$f"; fi;  done
 
#convert to a .nc4 file
for f in *.asc;

do
    year=`echo $f | cut -c16-19`
    month=`echo $f | cut -c20-21`
    day=`echo $f | cut -c22-23`
    outName=`echo $f | cut -c1-23`
    
    if [[ ! -f convert_2_nc/${outName}.nc4 ]]; then
       cdo -f nc -remapcon,$mask/CONUS_mask.grd -settaxis,${year}-${month}-${day},00:00:00,1day -setname,EDDI -input,$data_s/SMERGE_4_EDDI.grd  convert_2_nc/${outName}.nc4 < $f
       fi;
done

ulimit -n 8300
#Merge files

cdo mergetime convert_2_nc/*.nc4 convert_2_nc/eddi_merged.nc4



#### MERRA-2 data






