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
scripts=$data_s/NCAR_scripts
outData=$data_d/SubX/fromCasper
model_array=(GMAO ESRL RSMAS ECCC NRL EMC)


#Move data to HPC (easley)
rsync -Pa --ignore-existing  ~/Insync/OneDrive/NRT_CPC_Internship/ kdl0013@easley.auburn.edu:/home/kdl0013/NRT_CPC_Internship/
#Return data from HPC (easley)
rsync -PaL --ignore-existing kdl0013@easley.auburn.edu:/home/kdl0013/return_data_ETo/* ~/Insync/OneDrive/NRT_CPC_Internship/Data/SubX/
#random sending to easley
rsync -Pa --ignore-existing ~/Insync/OneDrive/NRT_CPC_Internship/Data/SubX/ECCC kdl0013@easley.auburn.edu:/home/kdl0013/NRT_CPC_Internship/Data/SubX/

###################################### STEP 1 ###################################
#Create wget scripts for each variable
#TODO:Return netcdf file forecasts from GEFSv12 from EASLEY####################
echo Mi
return_easley_files () {
#Return scripts for downloading and converting data
rsync -Pa kdl0013@easley.auburn.edu:/home/kdl0013/download_process_GEFSv12/GEFSv12/* ~/Insync/OneDrive/NRT_CPC_Internship/Data/SubX/*/raw_ensemble_files 
rsync -Pa kdl0013@easley.auburn.edu:/home/kdl0013/download_process_GEFSv12/process_GEFS/* ~/Insync/OneDrive/NRT_CPC_Internship/Scripts/EASLEY_HPC/ && \
rsync -Pa kdl0013@easley.auburn.edu:/home/kdl0013/download_process_GEFSv12/wget* ~/Insync/OneDrive/NRT_CPC_Internship/Scripts/EASLEY_HPC/
}
#Contains years 2020 - present from easley
return_easley_files #for GEFSv12 only

#maybe remove d35_dswrf_sfc_2018072500_p04.nc4          tmp_2m/d35_tmp_2m_2003031900_p08.nc4

ulimit -n 100000
echo An
#TODO:Return netcdf file forecasts from CASPER (all models, part of EMC)#######
return_CASPER_SUBx_files () {
mkdir -r $1
mkdir -r $2
rsync -Pa --ignore-existing klesinger@cheyenne.ucar.edu:/glade/u/home/klesinger/SubX/*/compress_resize/* $2
rsync -Pa --ignore-existing klesinger@cheyenne.ucar.edu:/glade/u/home/klesinger/SubX/S* $1
rsync -Pa --ignore-existing klesinger@cheyenne.ucar.edu:/glade/u/home/klesinger/SubX/wget* $1
}
return_CASPER_SUBx_files $scripts $outData
#This works better as of this point. Too many files
rsync -Pa --ignore-existing klesinger@cheyenne.ucar.edu:/glade/u/home/klesinger/SubX/ECCC/compress_resize/* /home/kdl/Insync/OneDrive/NRT_CPC_Internship/Data/SubX/fromCasper

#TODO: Remove s dimension for later processing
remove_s_dim () {
#ncks operators don't like overwriting in python without needing input from user; so don't compress
cat $data_s/01e_remove_S_dimsension_all_models.py | sed 's|model_name|'${1}'|g'> $data_so/$1_01e_remove_S_dimsension_all_models.py && python3 $data_so/$1_01e_remove_S_dimsension_all_models.py 
}
#model_array=(ESRL)
for model in "${model_array[@]}";
do remove_s_dim $model;done

python $data_s/fix_NRL_dates.py

#TODO: Preprocess files sent from ESRL FIMr1p1 agency #########################
{
ESRL_dir=$data_d/SubX/ESRL/from_ESRL_agency_missing_days
mkdir -p $ESRL_dir/keep_variables

cd $ESRL_dir

for zip in *.zip;do
unzip $zip;done

cp tas_2m* rad* keep_variables
rm *.nc #remove everything else for now

#restrict coordinates
cd $ESRL_dir/keep_variables

for f in *.nc;do 
cdo -f nc -remapcon,$mask/CONUS_mask.grd $f "$f"4;done

python3 $data_s/ESRL1_rename_files.py
}
###Preprocess data
#TODO: Make a script to only get the correct days for RSMAS#####################
python3 $data_s/01d_select_RSMAS_dates.py
#TODO: Create daily averages. Cannot do cdo operators because they split the files weirdly
python3 $data_s/EMC0_day_average_GEFSv12_HPC.py
python3 $data_s/EMC0_day_sum_GEFSv12_HPC.py
python3 $data_s/EMC1_merge_3_soil_moisture_fields_GEFSv12.py
#mv $data_d/SubX/EMC/raw_ensemble_files /media/kdl/Seagate_1/CPC_project/raw_ensemble_files #to save on storage space
#mv $data_d/SubX/EMC/raw_ensemble_files /media/kdl/Seagate_1/CPC_project/raw_ensemble_files2/*
########################## END OF STEP 1 ######################################

########################## STEP 2 #############################################

#####################




#Make ETo mean 
ETo_and_all_variables_anomaly_and_correlation () {
var_list=("tasmin" "tasmax" "tas" "actual_vapor_pressure" "srad" "windspeed" "pr")
#Create anomaly 
for var in "${var_list[@]}";do
cat $data_s/anomaly_compute_from_mean_all_variables.py | sed 's|model_name|'${1}'|g' | sed 's|var_name|'${var}'|g' > $data_so/$1_anomaly_compute_from_mean_$2.py && python3 $data_so/$1_anomaly_compute_from_mean_$2.py;done

}


model_array=(GMAO)
for model in "${model_array[@]}";
do ETo_and_all_variables_anomaly_and_correlation "$model" "Penman";
done

for var in "${var_list[@]}";do
cat $data_s/anomaly_correlation_all_other_variables.py | sed 's|model_name|'${1}'|g' | sed 's|var_name|'${var}'|g' > $data_so/$1_anomaly_correlation_"$var".py && python3 $data_so/$1_anomaly_correlation_"$var".py;done



model_array=(EMC)
#Make ETo mean 
ETo_and_all_variables_anomaly_and_correlation () {
cat $data_s/ETo_reference_ET_Penman.py | sed 's|model_name|'${1}'|g'> $data_so/"$1"_ETo_reference_ET_"{$2}".py && python3 $data_so/"$1"_ETo_reference_ET_"{$2}".py
#Make ETo anomaly mean files
cat $data_s/anomaly_mean_ETo_not_EMC.py | sed 's|model_name|'${1}'|g' | sed 's|evap_equation|'${2}'|g' > $data_so/$1_anomaly_mean_$2_not_EMC.py && python3 $data_so/$1_anomaly_mean_$2_not_EMC.py
cat $data_s/anomaly_mean_ETo_only_EMC_11_models.py | sed 's|model_name|'${1}'|g' | sed 's|evap_equation|'${2}'|g' > $data_so/$1_anomaly_mean_ETo_$2_11_models.py && python3 $data_so/$1_anomaly_mean_ETo_$2_11_models.py
#Make the observations in same format to go with anomaly correlation
cat $data_s/MERRA0_reformat_to_SubX_ETo.py | sed 's|model_name|'${1}'|g' | sed 's|evap_equation|'${2}'|g' > $data_so/$1_MERRA0_reformat_to_SubX_ETo_$2.py && python3 $data_so/$1_MERRA0_reformat_to_SubX_ETo_$2.py
#Create anomaly 
cat $data_s/anomaly_compute_from_mean_ETo.py | sed 's|model_name|'${1}'|g' | sed 's|evap_equation|'${2}'|g' > $data_so/$1_anomaly_compute_from_mean_ETo_$2.py && python3 $data_so/$1_anomaly_compute_from_mean_ETo_$2.py
#Plot correlation
cat $data_s/anomaly_correlation_ETo.py | sed 's|model_name|'${1}'|g' | sed 's|evap_equation|'${2}'|g' > $data_so/$1_anomaly_correlation_$2.py && python3 $data_so/$1_anomaly_correlation_$2.py

if [[  $1 != "NRL"  ]];then
cat $data_s/CRPS_climpred_skill_ETo.py | sed 's|model_name|'${1}'|g' | sed 's|evap_equation|'${2}'|g' > $data_so/$1_CRPS_climpred_skill_ETo_$2.py && python3 $data_so/$1_CRPS_climpred_skill_ETo_$2.py
fi

cat $data_s/ACC_climpred_skill_ETo.py | sed 's|model_name|'${1}'|g' | sed 's|evap_equation|'${2}'|g' > $data_so/$1_ACC_climpred_skill_ETo_$2.py && python3 $data_so/$1_ACC_climpred_skill_ETo_$2.py
}

model_array=(EMC ESRL ECCC GMAO ESRL RSMAS)
for model in "${model_array[@]}";
do ETo_and_all_variables_anomaly_and_correlation "$model" "Penman";
done


anomaly_all_other_variables() {
#Create individual variable anomalies and correlation
model_array=(ESRL GMAO RSMAS)

for model in "${model_array[@]}";do

var_list=("tas" "actual_vapor_pressure" "srad"  "windspeed" "tasmin" "tasmax" )

for var in "${var_list[@]}";do

cat $data_s/anomaly_mean_individual_variables.py | sed 's|model_name|'${model}'|g' | sed 's|var_name|'${var}'|g' > $data_so/"$model"_anomaly_mean_"$var".py && python3 $data_so/"$model"_anomaly_mean_"$var".py

#cat $data_s/anomaly_mean_individual_variables_single_day_not_EMC.py | sed 's|model_name|'${1}'|g' | sed 's|var_name|'${var}'|g' > $data_so/"$model"_anomaly_mean_"$var"_single_day_not_EMC.py && python3 $data_so/"$model"_anomaly_mean_"$var"_single_day_not_EMC.py

#cat $data_s/MERRA0_reformat_to_SubX_all_variables.py | sed 's|model_name|'${1}'|g' | sed 's|var_name|'${var}'|g' > $data_so/$1_MERRA0_reformat_to_SubX_"$var".py && python3 $data_so/$1_MERRA0_reformat_to_SubX_"$var".py

cat $data_s/anomaly_compute_from_mean_all_variables.py | sed 's|model_name|'${model}'|g' | sed 's|var_name|'${var}'|g' > $data_so/"$model"_anomaly_compute_from_mean_"$var".py && python3 $data_so/"$model"_anomaly_compute_from_mean_"$var".py

cat $data_s/ACC_climpred_skill_all_variables.py | sed 's|model_name|'${model}'|g' | sed 's|var_name|'${var}'|g' > $data_so/"$model"_ACC_climpred_skill_"$var".py && python3 $data_so/"$model"_ACC_climpred_skill_"$var".py;done;done

}
anomaly_all_other_variables
anomaly_all_other_variables


#FIRST SEND TO HPC, it is very costly because of the function (KEEP THIS PART)
make_ETo_Penman_single () {
cat $data_s/ETo_reference_ET_Penman.py | sed 's|model_name|'${1}'|g'> $data_so/"$1"_ETo_reference_ET_Penman.py && python3 $data_so/"$1"_ETo_reference_ET_Penman.py
}
make_ETo_Penman_single "ECCC"

#####################
#create the mean file (to subtract and make anomaly)
echo && echo && echo
echo "Should we remove ETo $2 ... (yes/no)?"
read
if [[ $REPLY == "yes" ]]
then
#rm $data_d/SubX/$1/anomaly/mean_for_ACC/*$2_*
fi


for model in "${model_array[@]}";
do  "$model";
done


make_anomaly_corr_plots() {
var_list=("tasmin" "tasmax" "actual_vapor_pressure" "srad" "windspeed" "tas")
for var in "${var_list[@]}";do
cat $data_s/ACC_climpred_skill_all_variables.py | sed 's|model_name|'${1}'|g' | sed 's|var_name|'${var}'|g' > $data_so/$1_ACC_climpred_skill_"$var".py && python3 $data_so/$1_ACC_climpred_skill_"$var".py;done

}

model_array=(GMAO ESRL RSMAS ECCC NRL)
for model in "${model_array[@]}";
do make_anomaly_corr_plots "$model";
done


#&& python3 $data_s/$1_MERRA0_reformat_to_SubX_"$var".py
#var_list=("tasmin" "tasmax" "tas" "actual_vapor_pressure" "srad" "windspeed" "pr")
var_list=("tasmin" "tasmax")
test "ECCC" "Penman"


'("tasmin" "tasmax" "tas" "actual_vapor_pressure" "srad" "windspeed" "pr")'



#All models in 1 script for ACC correlation. This is my own special script for 
#getting anomaly averages before the ACC score is computed
python3 $data_s/anomaly_correlation_ETo_all_models.py 



cat $data_s/MERRA0_reformat_to_SubX_all_variables.py | sed 's|var_list|'${3}'|g' | sed 's|model_name|'${}'|g' > $data_s/$




ETo_mean_anomaly_and_correlation "ESRL" "Priestley"
ETo_mean_anomaly_and_correlation "RSMAS" "Priestley"

ETo_mean_anomaly_and_correlation "GMAO" "Priestley"








EMC_anomaly_mean_RZSM () {
#Make ETo mean file

#Make RZSM mean file
cat $data_s/anomaly_mean_RZSM_only_EMC_11_models.py | sed 's|model_name|'${1}'|g' > $data_so/$1_anomaly_mean_RZSM_11_models.py && python3 $data_so/$1_anomaly_mean_RZSM_11_models.py

#RZSM
cat $data_s/MERRA0_reformat_to_SubX_RZSM.py | sed 's|model_name|'${1}'|g' > $data_so/$1_MERRA0_reformat_to_SubX_RZSM.py && python3 $data_so/$1_MERRA0_reformat_to_SubX_RZSM.py
cat $data_s/anomaly_compute_from_mean_RZSM.py | sed 's|model_name|'${1}'|g' > $data_so/$1_anomaly_compute_from_mean_RZSM.py && python3 $data_so/$1_anomaly_compute_from_mean_RZSM.py

cat $data_s/anomaly_correlation_ETo.py | sed 's|model_name|'${1}'|g' | sed 's|evap_equation|'${2}'|g' > $data_so/$1_anomaly_correlation_$2.py && python3 $data_so/$1_anomaly_correlation_$2.py

cat $data_s/anomaly_correlation_RZSM.py | sed 's|model_name|'${1}'|g'  > $data_so/$1_anomaly_correlation_RZSM.py && python3 $data_so/$1_anomaly_correlation_RZSM.py
}
EMC_anomaly_mean "EMC" "Penman"






EMC_anomaly_mean () {

#cat $data_s/ACC_skill_RZSM.py | sed 's|model_name|'${1}'|g' > $data_so/$1_ACC_skill_RZSM.py && python3 $data_so/$1_ACC_skill_RZSM.py
#cat $data_s/CRPS_skill_ETo.py | sed 's|model_name|'${1}'|g' | sed 's|evap_equation|'${2}'|g' > $data_so/$1_CRPS_skill_ETo_$2.py && python3 $data_so/$1_CRPS_skill_ETo_$2.py
cat $data_s/CRPS_skill_RZSM.py | sed 's|model_name|'${1}'|g' > $data_so/$1_CRPS_skill_RZSM.py && python3 $data_so/$1_CRPS_skill_RZSM.py
}

EMC_anomaly_mean "EMC" "Priestley"

EMC_anomaly_mean "EMC" "Penman"














































#TODO: Penman Monteith equation (make models individually at first)
python3 $data_s/ETo_reference_ET_GMAO_Penman.py
python3 $data_s/ETo_reference_ET_EMC_Penman.py


#Make ETo mean 
make_ETo_Penman_mean_for_anomaly () {
cat $data_s/anomaly_mean_Penman_not_EMC.py | sed 's|model_name|'${1}'|g'> $data_s/$1_anomaly_mean_Penman_not_EMC.py && python3 $data_s/$1_anomaly_mean_Penman_not_EMC.py
#Make the observations in same format to go with anomaly correlation
cat $data_s/MERRA0_reformat_to_SubX_Penman.py | sed 's|model_name|'${1}'|g'> $data_s/$1_MERRA0_reformat_to_SubX_Penman.py && python3 $data_s/$1_MERRA0_reformat_to_SubX_Penman.py
#Create anomaly 
#cat $data_s/anomaly_compute_from_mean_Penman.py | sed 's|model_name|'${1}'|g'> $data_s/$1_anomaly_compute_from_mean_Penman.py && python3 $data_s/$1_anomaly_compute_from_mean_Penman.py
#anomaly correlation
cat $data_s/anomaly_correlation_ETo.py | sed 's|model_name|'${1}'|g' |> $data_s/$1_anomaly_correlation_ETo.py && python3 $data_s/$1_anomaly_correlation_ETo.py
}
make_ETo_Penman_mean_for_anomaly "GMAO"




#Turn observations into SubX format

python3 $data_s/MERRA0_reformat_to_SubX_Priestley_RSMAS.py
python3 $data_s/MERRA0_reformat_to_SubX_Priestley_ESRL.py

#TODO: GMAO GEOS-5 process
python3 $data_s/ETo_reference_ET_GMAO_Priestley.py
python3 $data_s/ETo_reference_ET_GMAO_Penman.py
python3 $data_s/GMAO0_remove_unncessary_dates.py
python3 $data_s/GMAO1_make_anomaly_mean_Priestley.py #make anomaly and multi ensemble mean anomaly
python3 $data_s/GMAO1_make_anomaly_mean_Penman.py
python3 $data_s/MERRA0_reformat_to_SubX_Priestley_GMAO.py
python3 $data_s/MERRA0_reformat_to_SubX_Penman_GMAO.py
python3 $data_s/GMAO4_compute_anomaly_from_mean_Priestley.py
python3 $data_s/GMAO5_anomaly_correlation.py 


#TODO: ESRL FIMr1p1
#Process model data produced from ESRL FIMr1p1 seperatley (for me :) )
python3 $data_s/ESRL1_rename_files.py
python3 $data_s/ETo_reference_ET_ESRL_Priestley.py
python3 $data_s/ESRL2_make_anomaly_mean_Priestley.py
python3 $data_s/ESRL3_compute_anomaly_from_mean_Priestley.py


#TODO: RSMAS CCSM4
#Make evaporative demand datasets (All models)
python3 $data_s/ETo_reference_ET_RSMAS.py
python3 $data_s/RSMAS1_make_anomaly_mean_Priestley.py
python3 $data_s/RSMAS2_compute_anomaly_from_mean_Priestley.py
python3 $data_s/RSMAS3_anomaly_correlation_Priestley.py 



#TODO: EMC GEFSv12
#Get the daily average, cdo operators do not do this properly because because the files
#were originally split. First 10 lead days. The 10th day had its data split between 2 seperate files.
python3 $data_s/EMC0_day_average_GEFSv12_HPC.py
python3 $data_s/EMC1_merge_3_soil_moisture_fields_GEFSv12.py
python3 $data_s/ETo_reference_ET_EMC_Priestley.py
python3 $data_s/EMC2_make_anomaly_mean_Priestley.py

#TODO: need to make ETo for gefsv12


#Deangelis et al. (2020). Don't include these years from GEFSv12
rm $subx/EMC/*2011* $subx/EMC/*2012*




# EMC GEFSv12 Get the daily average, cdo operators do not do this properly because because the files
#were originally split. First 10 lead days. The 10th day had its data split between 2 seperate files.
python3 $data_s/01a_day_average_GEFSv12_HPC.py
python3 $data_s/01b_merge_3_soil_moisture_fields_GEFSv12.py


### At this point, all model data has the S dimension removed, and has 
#already been verified to have data. Also all years 2000-2022 are in the same directory



make_ETo_Priestley () {
python3 $data_s/ETo_reference_ET_GMAO_Priestley.py
python3 $data_s/ETo_reference_ET_ESRL_Priestley.py
python3 $data_s/ETo_reference_ET_RSMAS_Priestley.py
python3 $data_s/ETo_reference_ET_EMC_Priestley.py
}




#OLD STUFF
#TODO: Remove unncessary dates GMAO

