#!/bin/bash

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

#Create wget scripts for each variable
#TODO:Return netcdf file forecasts from GEFSv12 from EASLEY
echo Mi
return_easley_files () {
rsync -Pa kdl0013@easley.auburn.edu:/home/kdl0013/GEFSv12/* ~/Insync/OneDrive/NRT_CPC_Internship/Data/SubX/EMC/raw_ensemble_files/ 

#Return scripts for downloading and converting data
rsync -Pa kdl0013@easley.auburn.edu:/home/kdl0013/process_GEFS/* ~/Insync/OneDrive/NRT_CPC_Internship/Scripts/EASLEY_HPC/ && \
rsync -Pa kdl0013@easley.auburn.edu:/home/kdl0013/wget* ~/Insync/OneDrive/NRT_CPC_Internship/Scripts/EASLEY_HPC/
}

return_easley_files

# EMC GEFSv12 Get the daily average, cdo operators do not do this properly because because the files
#were originally split. First 10 lead days. The 10th day had its data split between 2 seperate files.
python3 $data_s/01a_day_average_GEFSv12_HPC.py
python3 $data_s/01b_merge_3_soil_moisture_fields_GEFSv12.py


echo An
#TODO:Return netcdf file forecasts from CASPER (all models, part of EMC)
return_CASPER_SUBx_files () {
mkdir -r $1
mkdir -r $2
rsync -Pa klesinger@cheyenne.ucar.edu:/glade/u/home/klesinger/SubX/S* $1
rsync -Pa klesinger@cheyenne.ucar.edu:/glade/u/home/klesinger/SubX/*/compress_resize/* $2
rsync -Pa klesinger@cheyenne.ucar.edu:/glade/u/home/klesinger/SubX/wget* $1
}
return_CASPER_SUBx_files $scripts $outData
#Return file format ----- EMC.nc, GMAO.nc,



###Preprocess data
#TODO: Make a script to only get the correct days for RSMAS -- DONE
python3 $data_s/01d_select_RSMAS_dates.py

#TODO: Remove s dimension for later processing
python3 $data_s/01e_remove_S_dimsension_all_models.py

#TODO: Move files into seperate directories
model_array=(GMAO ESRL RSMAS EMC)

for model in "${model_array[@]}";
do mkdir $subx/$model
cp $outData/S_dim_removed/*$model*.nc4 $subx/$model/;
done



### At this point, all model data has the S dimension removed, and has 
#already been verified to have data. Also all years 2000-2022 are in the same directory

