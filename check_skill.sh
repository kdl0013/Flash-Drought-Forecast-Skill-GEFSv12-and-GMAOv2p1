#!/bin/bash

main_directory='/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
processors=8

data_d=$main_directory/Data
data_s=$main_directory/Scripts
data_so=$data_s/model_scripts
subx=$data_d/SubX
mask=$data_s/CONUS_mask


#SubX GEOS-5, FIMr1p1, CCSM4 hindcasts/forecasts, and ESRL forecasts
#return from Cheyenne UCAR cluster
outData=$data_d/SubX
model_array=(GMAO ESRL RSMAS ECCC NRL EMC)

for mod in "${model_array[@]}";do
echo $mod && sleep 1
ls $outData/$mod/anomaly/skill_assessments/*
echo
ls $outData/$mod/anomaly/skill_assessments/* | wc -l 
echo total vars completed out of 8 maximum
echo "##########################################################"

echo
echo 'All vars "ETo" "tasmax" "actual_vapor_pressure" "srad" "windspeed" "tasmin" "tas"'

sleep 2;done
