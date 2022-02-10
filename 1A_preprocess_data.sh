#!/bin/bash

#Download SMERGE soil moisture, gridMET PET, EDDI evaporative demand, and elevation
#To download SubX data, use script download_SubX_4_models.sh (run in HPC)

module load cdo
module load nco

#MUST SET inputs: main_directory and number of processors for CPU parallelism

#You can make main_directory path any path you would like and this will get
#all other directories in order. 
main_directory='/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
processors=8

data_d=$main_directory/Data
data_s=$main_directory/Scripts

############## Reference ETo data creation ###########
cat $data_s/1a_GMAO_compute_ET0_single_file.py | sed 's|main_dir|'${main_directory}'|g' | \
    sed 's|procs|'${processors}'|g' > $data_s/1a_TMP_GMAO_compute_ET0_single_file.py

echo Starting on daily reference ET calculation

python3 $data_s/1a_TMP_GMAO_compute_ET0_single_file.py
rm $data_s/1a_TMP_GMAO_compute_ET0_single_file.py

############## Reference ETo and Soil Moisture scatterplot data creation ###########
cat $data_s/1b_GMAO_scatterplots_ET0_RZSM.py | sed 's|main_dir|'${main_directory}'|g' | \
    sed 's|procs|'${processors}'|g' > $data_s/1b_TMP_GMAO_scatterplots_ET0_RZSM.py

echo Starting on daily reference ET calculation

python3 $data_s/1b_TMP_GMAO_scatterplots_ET0_RZSM.py
rm $data_s/1b_TMP_GMAO_scatterplots_ET0_RZSM.py

############## EDDI pre-processing ###########








############## Lagged Average Ensemble ##############
cat $data_s/3a_lagged_average_ensemble.py | sed 's|main_dir|'${main_directory}'|g' > $data_s/3a_TMP_lagged_average_ensemble.py

python3 $data_s/3a_TMP_lagged_average_ensemble.py
rm $data_s/3a_TMP_lagged_average_ensemble.py










