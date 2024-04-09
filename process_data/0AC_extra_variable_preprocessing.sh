#!/bin/bash
cas #load spyder and other packages
main_directory='/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
processors=8

data_d=$main_directory/Data
data_s=$main_directory/Scripts
subx=$data_d/SubX
mask=$data_s/CONUS_mask


#SubX GEOS-5, FIMr1p1, CCSM4 hindcasts/forecasts, and ESRL forecasts
#return from Cheyenne UCAR cluster
scripts=$data_s/NCAR_scripts
outData=$data_d/SubX/fromCasper
model_array=(GMAO RSMAS EMC ECCC NRL ESRL)

#Check latest dates downloaded
for mod in "${model_array[@]}";do
ls -ls $data_d/SubX/$mod | tail -7;done

#Priestley-Taylor coefficient -- COMPLETE
pt_dir=$data_d/Priestley_Taylor_Makkink_evap_coeff_maps
ncatted -O -a units,X,c,c,"degrees_east" -a units,Y,c,c,"degrees_south" pt_coeff.nc pt_coeff_add_metadata.nc

cdo -remapbil,$mask/CONUS_mask.grd $pt_dir/pt_coeff_add_metadata.nc $pt_dir/pt_coeff_FINAL.nc
################################################################################
#Elevation dataset preprocess (for altitude) -- COMPLETE
mkdir $data_d/elevation && cd $data_d/elevation
data_e=$data_d/elevation

cdo -remapbil,$mask/CONUS_mask.grd $data_e/elev.1-deg.nc $data_e/elev_regrid.nc
################################################################################


mkdir $data_d/elevation && cd $data_d/elevation
data_e=$data_d/elevation

cdo -remapbil,$mask/CONUS_mask.grd $data_e/elev.1-deg.nc $data_e/elev_regridtest.nc







