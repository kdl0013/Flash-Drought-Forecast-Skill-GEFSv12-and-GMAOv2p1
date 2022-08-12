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
model=GEFSv12
processors=7

data_d=$main_directory/Data
data_s=$main_directory/Scripts
sub_dir=$data_d/SubX/$model
mkdir -p $sub_dir
mkdir -p $data_s 
mkdir -p $data_d/GLEAM_RZSM
mkdir -p $data_d/gridMET/ETo_SubX_values
#mkdir -p $data_d/EDDI/convert_2_nc
################### SMERGE soil moisture ##########################

cd $data_s
#Create wget scripts for each variable
#Return netcdf file forecasts from GEFSv12 from EASLEY
return_easley_files () {
rsync -Pa kdl0013@easley.auburn.edu:/home/kdl0013/GEFSv12/* ~/Insync/OneDrive/NRT_CPC_Internship/Data/SubX/GEFSv12/raw_ensemble_files/ 

#Return scripts for downloading and converting data
rsync -Pa kdl0013@easley.auburn.edu:/home/kdl0013/process_GEFS/multi* ~/Insync/OneDrive/NRT_CPC_Internship/Scripts/EASLEY_HPC/ && \
rsync -Pa kdl0013@easley.auburn.edu:/home/kdl0013/wget* ~/Insync/OneDrive/NRT_CPC_Internship/Scripts/EASLEY_HPC/
}

return_easley_files


ensem_dirs=$sub_dir/raw_ensemble_files
cd $ensem_dirs
#TODO: Since the files are already day average we need to combine each model into
#one file (like SubX format)

python3 01_combine_models_GEFSv12.py
    









