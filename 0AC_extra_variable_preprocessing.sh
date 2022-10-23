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

cdo -remapcon,$mask/CONUS_mask.grd $pt_dir/pt_coeff_add_metadata.nc $pt_dir/pt_coeff_FINAL.nc
################################################################################
#Elevation dataset preprocess (for altitude) -- COMPLETE
mkdir $data_d/elevation && cd $data_d/elevation
data_e=$data_d/elevation

cdo -remapcon,$mask/CONUS_mask.grd $data_e/elev.1-deg.nc $data_e/elev_regrid.nc
################################################################################


#TODO: Test GMAO radiation
python3 $data_s/MERRA0_reformat_to_SubX_radiation.py
python3 $data_s/anomaly_correlation_radiation.py


#TODO: Create individual variable anomalies
anomaly_by_variable () {
cat $data_s/anomaly_mean_individual_variables_not_EMC.py | sed 's|model_name|'${1}'|g' | sed 's|var_name|'${2}'|g' > $data_s/$1_anomaly_mean_$2_not_EMC.py && python3 $data_s/$1_anomaly_mean_$2_not_EMC.py

}


vars=("pr" "tasmin" "tasmax" "tas" "actual_vapor_pressure" "srad" "windspeed")
for model in "${model_array[@]}";do
for var in "${vars[@]}";do 
anomaly_by_variable "$model" "$var";done;done


RSMAS_vars=("tasmin" "tasmax" "tas" "srad" "windspeed" "pr" "actual_vapor_pressure")
for var in "${RSMAS_vars[@]}";
do anomaly_by_variable "RSMAS" "$var";done









anomaly_by_variable "GMAO" "tasmax"
anomaly_by_variable "GMAO" "dswrf"
anomaly_by_variable "GMAO" "dlwrf"
anomaly_by_variable "GMAO" "uswrf"
anomaly_by_variable "GMAO" "ulwrf"
anomaly_by_variable "GMAO" "tdps"
anomaly_by_variable "GMAO" "uas"
anomaly_by_variable "GMAO" "vas"
anomaly_by_variable "GMAO" "srad"



####### MERRA2 make anomalies
python3 $data_s/OBS_make_anomalies_Penman.py
python3 $data_s/OBS_make_anomalies_Priestley_and_RZSM.py



























#Find the mean for each file first
python3 $data_s/ALL_MODEL_EMC_make_anomaly_mean.py

python3 $data_s/ALL_MODEL_GMAO_make_anomaly_mean.py

python3 $data_s/ALL_MODEL_RSMAS_make_anomaly_mean.py

python3 $data_s/ALL_MODEL_ESRL_make_anomaly_mean.py

#Because of some memory issues, I need to process some of them in the cluster
rsync -Pa $data_d/SubX/* klesinger@cheyenne.ucar.edu:/glade/scratch/klesinger/SubX_process
#scripts only
rsync -Pa $data_s/* klesinger@cheyenne.ucar.edu:/glade/u/home/klesinger/scripts_linux
#CONUs mask only
rsync -Pa $data_d/CONUS_mask/* klesinger@cheyenne.ucar.edu:/glade/u/home/klesinger/scripts_linux

#OLD
#GMAO Soil moisture. Convert depth file given by Dr. Robert Koster to netcdf (python)
python3 $data_s/GMAO0_convert_depth_dat_file.py



#Convert soil moisture GMAO to proper format (m3/m3)
python3 $data_s/GMAO2_convert_RZSM_to_m3_m3.py

#Make evaporative demand dataset GMAO
python3 $data_s/GMAO3_make_reference_ET.py

#Create a mean dataset for anomaly
python3 $data_s/GMAO4_make_anomaly_mean.py


