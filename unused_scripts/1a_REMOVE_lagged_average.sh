#!/bin/bash

#Download SMERGE soil moisture, gridMET PET, and EDDI evaporative demand
#To download SubX data, use script download_SubX_4_models.sh (run in HPC)

module load cdo
module load nco

#MUST SET inputs: main_directory and number of processors for CPU parallelism

###SMERGE soil moisture
#You can make main_directory path any path you would like. 
#This path will allow everything to download in proper order
main_directory='/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
processors=12

data_d=$main_directory/Data
data_s=$main_directory/Scripts
mkdir -p $data_d/Data/SubX/GMAO && mkdir -p $data_s

##### Lagged Average Ensemble #####
cat 1a_lagged_average_ensemble.py | sed 's|main_dir|'${main_directory}'|g' > 1a_TMP_lagged_average_ensemble.py

python3 $data_s/1a_TMP_lagged_average_ensemble.py && rm 1a_TMP_lagged_average_ensemble.py


