#!/bin/bash

source ~/.bashrc
cas


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
var_list_ETo=("ETo_Penman" "tas" "actual_vapor_pressure" "total_radiation" "windspeed" "pr" "zg")
var_list=("tasmin" "tasmax" "tas" "actual_vapor_pressure" "total_radiation" "windspeed" "pr" "zg")


########################## STEP 2 #############################################


#Make ETo mean, compute anomaly, find correlation between reanalysis 
ETo_and_all_variables_anomaly_and_correlation () {
cat $data_s/ETo_reference_ET_Penman_FAO56_reforecast_EMC.py 

#Make ETo anomaly mean files (using only shortwave radiation)
#cat $data_s/anomaly_mean_ETo_not_EMC_shortwave.py | sed 's|model_name|'${1}'|g' | sed 's|evap_equation|'${2}'|g' > $data_so/$1_anomaly_mean_$2_not_EMC.py && python3 $data_so/$1_anomaly_mean_$2_not_EMC.py

#Make ETo anomaly mean files (using net radiation calculated from extraterrestrial, and proper reference ET conditions
cat $data_s/anomaly_mean_ETo_not_EMC_total_radiation.py | sed 's|model_name|'${1}'|g' | sed 's|evap_equation|'${2}'|g' > $data_so/$1_anomaly_mean_$2_not_EMC.py && python3 $data_so/$1_anomaly_mean_$2_not_EMC.py

#cat $data_s/anomaly_mean_ETo_only_EMC_11_models.py | sed 's|model_name|'${1}'|g' | sed 's|evap_equation|'${2}'|g' > $data_so/$1_anomaly_mean_ETo_$2_11_models.py && python3 $data_so/$1_anomaly_mean_ETo_$2_11_models.py

#Make the observations in same format to go with anomaly correlation
cat $data_s/MERRA0_reformat_to_SubX_ETo.py | sed 's|model_name|'${1}'|g' | sed 's|evap_equation|'${2}'|g' > $data_so/$1_MERRA0_reformat_to_SubX_ETo_$2.py && python3 $data_so/$1_MERRA0_reformat_to_SubX_ETo_$2.py

#Make the observations in same format to go with anomaly correlation
cat $data_s/GEFS0_reformat_to_SubX_ETo.py | sed 's|model_name|'${1}'|g' | sed 's|evap_equation|'${2}'|g' > $data_so/$1_GEFS0_reformat_to_SubX_ETo_$2.py && python3 $data_so/$1_GEFS0_reformat_to_SubX_ETo_$2.py

#Make the observations in same format to go with anomaly correlation
cat $data_s/NLDAS0_reformat_to_SubX_ETo.py | sed 's|model_name|'${1}'|g' | sed 's|evap_equation|'${2}'|g' > $data_so/$1_GEFS0_reformat_to_SubX_ETo_$2.py && python3 $data_so/$1_GEFS0_reformat_to_SubX_ETo_$2.py

#Create anomaly 
cat $data_s/anomaly_compute_from_mean_ETo.py | sed 's|model_name|'${1}'|g' | sed 's|evap_equation|'${2}'|g' > $data_so/$1_anomaly_compute_from_mean_ETo_$2.py && python3 $data_so/$1_anomaly_compute_from_mean_ETo_$2.py

#Plot correlation
cat $data_s/anomaly_correlation_ETo_MERRA_GEFS.py | sed 's|model_name|'${1}'|g' | sed 's|evap_equation|'${2}'|g' | sed 's|reanalysis_select|MERRA|g' > $data_so/$1_anomaly_correlation_$2_MERRA.py && python3 $data_so/$1_anomaly_correlation_$2_MERRA.py

cat $data_s/anomaly_correlation_ETo_MERRA_GEFS.py | sed 's|model_name|'${1}'|g' | sed 's|evap_equation|'${2}'|g' | sed 's|reanalysis_select|GEFSv12|g' > $data_so/$1_anomaly_correlation_$2_GEFS.py && python3 $data_so/$1_anomaly_correlation_$2_GEFS.py

cat $data_s/anomaly_correlation_ETo_MERRA_GEFS.py | sed 's|model_name|'${1}'|g' | sed 's|evap_equation|'${2}'|g' | sed 's|reanalysis_select|NLDAS|g' > $data_so/$1_anomaly_correlation_$2_GEFS.py && python3 $data_so/$1_anomaly_correlation_$2_GEFS.py
#if [[  $1 != "NRL"  ]];then
#cat $data_s/CRPS_climpred_skill_ETo.py | sed 's|model_name|'${1}'|g' | sed 's|evap_equation|'${2}'|g' > $data_so/$1_CRPS_climpred_skill_ETo_$2.py && python3 $data_so/$1_CRPS_climpred_skill_ETo_$2.py
#fi

#Can verify data with climpred
#cat $data_s/ACC_climpred_skill_ETo.py | sed 's|model_name|'${1}'|g' | sed 's|evap_equation|'${2}'|g' > $data_so/$1_ACC_climpred_skill_ETo_$2.py && python3 $data_so/$1_ACC_climpred_skill_ETo_$2.py
}

#model_array=(ECCC NRL ESRL EMC GMAO RSMAS)
model_array=(EMC GMAO) #dont do any other models. They dont' have net radiation partitioned

#Only look at EMC and GMAO because they have downwelling shortwave radiation partitioned
for model in "${model_array[@]}";
do ETo_and_all_variables_anomaly_and_correlation "$model" "Penman";
done


python $data_s/check_bad_anomaly_mean.py #for anomaly mean files with bad outputs

anomaly_all_other_variables() {
#DONT DO EMC zg, there isn't enough data for a mean
#Create individual variable anomalies and correlation
model_array=(EMC ECCC ESRL RSMAS GMAO NRL)
#model_array=(EMC)
#model_array=(NRL)
for model in "${model_array[@]}";do

var_list=("tasmin" "tasmax" "tas" "actual_vapor_pressure" "total_radiation" "windspeed" "pr" "zg200" "dswrf")
var_list=('tasmin' 'tasmax')
for var in "${var_list[@]}";do

cat $data_s/anomaly_mean_individual_variables.py | sed 's|model_name|'${model}'|g' | sed 's|var_name|'${var}'|g' > $data_so/"$model"_anomaly_mean_"$var".py && python3 $data_so/"$model"_anomaly_mean_"$var".py

cat $data_s/MERRA0_reformat_to_SubX_all_variables.py | sed 's|model_name|'${model}'|g' | sed 's|var_name|'${var}'|g' > $data_so/"$model"_MERRA0_reformat_to_SubX_"$var".py && python3 $data_so/"$model"_MERRA0_reformat_to_SubX_"$var".py

cat $data_s/GEFS0_reformat_to_SubX_all_variables.py | sed 's|model_name|'${model}'|g' | sed 's|var_name|'${var}'|g' > $data_so/"$model"_GEFS0_reformat_to_SubX_"$var".py && python3 $data_so/"$model"_GEFS0_reformat_to_SubX_"$var".py

cat $data_s/NLDAS0_reformat_to_SubX_all_variables.py | sed 's|model_name|'${model}'|g' | sed 's|var_name|'${var}'|g' > $data_so/"$model"_GEFS0_reformat_to_SubX_"$var".py && python3 $data_so/"$model"_GEFS0_reformat_to_SubX_"$var".py

cat $data_s/anomaly_compute_from_mean_all_variables.py | sed 's|model_name|'${model}'|g' | sed 's|var_name|'${var}'|g' > $data_so/"$model"_anomaly_compute_from_mean_"$var".py && python3 $data_so/"$model"_anomaly_compute_from_mean_"$var".py

#with MERRA only
cat $data_s/anomaly_correlation_all_variables_MERRA_GEFS.py | sed 's|model_name|'${model}'|g' | sed 's|var_name|'${var}'|g' | sed 's|reanalysis_select|MERRA|g'> $data_so/"$model"_anomaly_correlation_"$var"_MERRA.py && python3 $data_so/"$model"_anomaly_correlation_"$var"_MERRA.py

#With GEFSv12 only
cat $data_s/anomaly_correlation_all_variables_MERRA_GEFS.py | sed 's|model_name|'${model}'|g' | sed 's|var_name|'${var}'|g' | sed 's|reanalysis_select|GEFSv12|g'> $data_so/"$model"_anomaly_correlation_"$var"_GEFS.py && python3 $data_so/"$model"_anomaly_correlation_"$var"_GEFS.py;done;done

#With NLDAS only
cat $data_s/anomaly_correlation_all_variables_MERRA_GEFS.py | sed 's|model_name|'${model}'|g' | sed 's|var_name|'${var}'|g' | sed 's|reanalysis_select|NLDAS|g'> $data_so/"$model"_anomaly_correlation_"$var"_GEFS.py && python3 $data_so/"$model"_anomaly_correlation_"$var"_GEFS.py;done;done

#This is only to see if my methodology works and it did for MERRA.
#cat $data_s/ACC_climpred_skill_all_variables_MERRA.py | sed 's|model_name|'${model}'|g' | sed 's|var_name|'${var}'|g' | sed 's|reanalysis_select|MERRA|g' > $data_so/"$model"_ACC_climpred_skill_"$var".py && python3 $data_so/"$model"_ACC_climpred_skill_"$var".py

}
anomaly_all_other_variables

#This is only to see if my methodology using just the seasonal mean (add arguments: Subx model, variable, reanlaysis model) with CLIMPRED package
model_array=("EMC" "GMAO")
var_list=("actual_vapor_pressure" "total_radiation")
reanalysis_list=("MERRA" "GEFSv12", "NLDAS")

for model in "${model_array[@]}";do
for var in "${var_list[@]}";do
for reanalysis in "${reanalysis_list[@]}";do

python3 $data_s/ACC_climpred_skill_all_variables_MERRA_GEFS.py $model $var $reanalysis;done;done;done


var_list=("total_radiation")
for model in "${model_array[@]}";do
for var in "${var_list[@]}";do
for reanalysis in "${reanalysis_list[@]}";do

python3 $data_s/ACC_climpred_skill_all_variables_normal_values_MERRA_GEFS.py $model $var $reanalysis;done;done;done

python3 $data_s/ACC_climpred_skill_all_variables_normal_values_MERRA_GEFS_actual_vapor_pressure_only.py

#TODO:CPC precipitation skill check (not doing this anymore)
#{
#model_array=(GMAO NRL ESRL ECCC RSMAS EMC)
#model_array=(ECCC ESRL)
#var='pr'
#for model in "${model_array[@]}";do
#cat $data_s/CPC_reformat_to_SubX_all_variables.py | sed 's|model_name|'${model}'|g' | sed 's|var_name|'${var}'|g' > #$data_so/"$model"_CPC_reformat_to_SubX_"$var".py && python3 $data_so/"$model"_CPC_reformat_to_SubX_"$var".py
#cat $data_s/anomaly_correlation_CPC_obs.py | sed 's|model_name|'${model}'|g' | sed 's|var_name|'${var}'|g' > $data_so/"$model"_anomaly_correlation_"$var"_CPC.py && python3 $data_so/"$model"_anomaly_correlation_"$var"_CPC.py;done
#}


#########All models have been processed, now combine them to make a single plot for each reanalysis
var_list=("ETo_Penman" "tas" "actual_vapor_pressure" "total_radiation" "pr" "windspeed" "zg200" "tasmax" "tasmin" "dswrf")

#var_list=("tas" "tasmin" "tasmax")
#for var in "${var_list[@]}";do

#cat $data_s/plot_correlation_ALL_MODEL_MERRA_GEFSv12.py | sed 's|var_name|'${var}'|g' | sed 's|reanalysis_select|GEFSv12|g' > $data_so/"$var"_ACC_correlation_ALL_MODEL_GEFSv12.py && python3 $data_so/"$var"_ACC_correlation_ALL_MODEL_GEFSv12.py

#cat $data_s/plot_correlation_ALL_MODEL_MERRA_GEFSv12.py | sed 's|var_name|'${var}'|g' | sed 's|reanalysis_select|MERRA|g' > $data_so/"$var"_ACC_correlation_ALL_MODEL_MERRA.py && python3 $data_so/"$var"_ACC_correlation_ALL_MODEL_MERRA.py;done

#cpc precipitation
#python3 $data_s/plot_correlation_ALL_MODEL_cpc.py

#var_list=("tas" "tasmin" "tasmax")
for var in "${var_list[@]}";do
cat $data_s/plot_correlation_ALL_MODEL_MERRA_GEFSv12_both_models.py | sed 's|var_name|'${var}'|g' > $data_so/"$var"_ACC_correlation_ALL_MODEL_GEFSv12_both_models.py && python3 $data_so/"$var"_ACC_correlation_ALL_MODEL_GEFSv12_both_models.py;done


EMC_anomaly_mean_RZSM () {
#Make ETo mean file

#Make RZSM mean file
cat $data_s/anomaly_mean_RZSM_only_EMC_11_models.py | sed 's|model_name|'${1}'|g' > $data_so/$1_anomaly_mean_RZSM_11_models.py && python3 $data_so/$1_anomaly_mean_RZSM_11_models.py

cat $data_s/anomaly_mean_RZSM_weighted_only_EMC_11_models.py | sed 's|model_name|'${1}'|g' > $data_so/$1_anomaly_mean_RZSM_weighted_11_models.py && python3 $data_so/$1_anomaly_mean_RZSM_weighted_11_models.py

cat $data_s/anomaly_compute_from_mean_RZSM.py | sed 's|model_name|'${1}'|g' > $data_so/$1_anomaly_compute_from_mean_RZSM.py && python3 $data_so/$1_anomaly_compute_from_mean_RZSM.py  #Only look at the weighted one, I've already run the other ones

#convert RZSM reanalysis and RZSM percentiles to same as SubX
cat $data_s/MERRA0_reformat_to_SubX_RZSM.py | sed 's|model_name|'${1}'|g' > $data_so/$1_MERRA0_reformat_to_SubX_RZSM.py && python3 $data_so/$1_MERRA0_reformat_to_SubX_RZSM.py

cat $data_s/GEFS_reformat_to_SubX_RZSM_weighted.py | sed 's|model_name|'${1}'|g' > $data_so/$1_GEFS_reformat_to_SubX_RZSM.py && python3 $data_so/$1_GEFS_reformat_to_SubX_RZSM.py

cat $data_s/NLDAS_reformat_to_SubX_RZSM.py | sed 's|model_name|'${1}'|g' > $data_so/$1_GEFS_reformat_to_SubX_RZSM.py && python3 $data_so/$1_GEFS_reformat_to_SubX_RZSM.py

#find smpd index on observations

python3 $data_s/1n_RZSM_FD_weighted_classification_reformatted_GEFSV12_renalysis.py #done


#create the SMPD index on EMC 
#added a 7-day running mean to perceniles first
python3 $data_s/1n_RZSM_FD_classification_SubX_all_models_EMC.py #done
python3 $data_s/1n_RZSM_FD_weighted_classification_SubX_all_models_EMC.py #done
#python3 $data_s/1n_RZSM_FD_classification_SubX_ensemble_mean_EMC.py #OLDj aoj




#cat $data_s/GEFS_reformat_to_SubX_RZSM_actual_values_no_anomalies.py | sed 's|model_name|'${1}'|g' > $data_so/$1_GEFS_reformat_to_SubX_RZSM_actual_values_no_anomalies.py && python3 $data_so/$1_GEFS_reformat_to_SubX_RZSM_actual_values_no_anomalies.py

#cat $data_s/anomaly_correlation_RZSM_MERRA_warm_season.py > $data_so/$1_anomaly_correlation_RZSM_MERRA_warm_season.py && python3 $data_so/$1_anomaly_correlation_RZSM_MERRA_warm_season.py

#cat $data_s/anomaly_correlation_RZSM_GEFS_warm_season.py > $data_so/$1_anomaly_correlation_RZSM_GEFS_warm_season.py && python3 $data_so/$1_anomaly_correlation_RZSM_GEFS_warm_season.py

cat $data_s/anomaly_correlation_RZSM_MERRA_all_season.py > $data_so/$1_anomaly_correlation_RZSM_MERRA_all_season.py && python3 $data_so/$1_anomaly_correlation_RZSM_MERRA_all_season.py

cat $data_s/anomaly_correlation_RZSM_GEFS_all_season.py > $data_so/$1_anomaly_correlation_RZSM_GEFS_all_season.py && python3 $data_so/$1_anomaly_correlation_RZSM_GEFS_all_season.py

cat $data_s/anomaly_correlation_RZSM_NLDAS_all_season.py > $data_so/$1_anomaly_correlation_RZSM_GEFS_all_season.py && python3 $data_so/$1_anomaly_correlation_RZSM_GEFS_all_season.py

cat $data_s/anomaly_correlation_RZSM_weighted_MERRA_all_season.py > $data_so/$1_anomaly_correlation_RZSM_weighted_MERRA_all_season.py && python3 $data_so/$1_anomaly_correlation_RZSM_weighted_MERRA_all_season.py

cat $data_s/anomaly_correlation_RZSM_weighted_GEFS_all_season.py > $data_so/$1_anomaly_correlation_RZSM_weighted_GEFS_all_season.py && python3 $data_so/$1_anomaly_correlation_RZSM_weighted_GEFS_all_season.py

cat $data_s/anomaly_correlation_RZSM_weighted_NLDAS_all_season.py > $data_so/$1_anomaly_correlation_RZSM_weighted_MERRA_all_season.py && python3 $data_so/$1_anomaly_correlation_RZSM_weighted_MERRA_all_season.py

cat $data_s/anomaly_correlation_RZSM_weighted_GEFS_all_gridcells.py > $data_so/$1_anomaly_correlation_RZSM_weighted_GEFS_all_gridcells.py && python3 $data_so/$1_anomaly_correlation_RZSM_weighted_GEFS_all_gridcells.py #don't need to do MERRA

cat $data_s/anomaly_correlation_RZSM_GEFS_all_gridcells.py > $data_so/$1_anomaly_correlation_RZSM_GEFS_all_gridcells.py && python3 $data_so/$1_anomaly_correlation_RZSM_GEFS_all_gridcells.py #don't need to do MERRA



#cat $data_s/RZSM_GEFS_day1_bias.py > $data_so/$1_RZSM_GEFS_day1_bias.py && python3 $data_so/$1_RZSM_GEFS_day1_bias.py

cat $data_s/RZSM_GEFS_day1_anomaly_bias.py > $data_so/$1_RZSM_GEFS_day1_anomaly_bias.py && python3 $data_so/$1_RZSM_GEFS_day1_anomaly_bias.py

cat $data_s/RZSM_weighted_GEFS_day1_anomaly_bias.py > $data_so/$1_RZSM_weighted_GEFS_day1_anomaly_bias.py && python3 $data_so/$1_RZSM_weighted_GEFS_day1_anomaly_bias.py #done

#check RZSM skill with predicted vs. actual precipitation
#cat $data_s/anomaly_correlation_RZSM_check_skill_of_dry_anomalies.py | sed 's|model_name|'${1}'|g'  > $data_so/$1_anomaly_correlation_RZSM_wet_dry_anomaly.py && python3 $data_so/$1_anomaly_correlation_RZSM_wet_dry_anomaly.py

#check RZSM skill by different categories (dryer vs. wetter RZSM anomalies)
#cat $data_s/anomaly_correlation_RZSM_check_skill_greater_less.py | sed 's|model_name|'${1}'|g'  > $data_so/$1_anomaly_correlation_RZSM_check_skill_greater_less.py && python3 $data_so/$1_anomaly_correlation_RZSM_check_skill_greater_less.py

cat $data_s/anomaly_correlation_RZSM_weighted_check_skill_greater_less.py | sed 's|model_name|'${1}'|g'  > $data_so/$1_anomaly_correlation_RZSM_weighted_check_skill_greater_less.py && python3 $data_so/$1_anomaly_correlation_RZSM_weighted_check_skill_greater_less.py

cat $data_s/plot_correlation_RZSM_all_wet_dry.py | sed 's|model_name|'${1}'|g' > $data_so/RZSM_ACC_correlation.py && python3 $data_so/RZSM_ACC_correlation.py;done

python3 $data_s/plot_correlation_RZSM_MERRA_GEFSv12.py

}
EMC_anomaly_mean "EMC" "Penman"



### check why radiation skill is different in the fall
python3 $data_s/anomaly_correlation_individual_variables_GEFS_all_gridcells_RZSM.py #run for RZSM and ETo_Penman



reanalysis_list=("MERRA" "GEFSv12", "NLDAS")

### Create the SMPD index on each variable
for reanalysis in "${reanalysis_list[@]}";do
python3 $data_s/1n_RZSM_FD_classification_reformatted_NLDAS_renalysis.py $reanalysis
python3 $data_s/SMPD_heidke_skill.py $reanalysis
python3 $data_s/SMPD_weighted_heidke_skill_MEM.py $reanalysis
python3 $data_s/plot_correlation_heidke_hitRate_weighted.py $reanalysis
python3 $data_s/plot_correlation_heidke_hitRate.py $reanalysis;
done

# For individual realizations (see code below)

#cat $data_s/SMPD_weighted_heidke_skill_individual_realization.py | sed 's|reanalysis_select|GEFSv12|g' > $data_so/SMPD_weighted_heidke_skill_individual_realization.py && python3 $data_so/SMPD_weighted_heidke_skill_individual_realization.py

#cat $data_s/SMPD_non_weighted_heidke_skill_individual_realization.py | sed 's|reanalysis_select|GEFSv12|g' > $data_so/SMPD_non_weighted_heidke_skill_individual_realization.py && python3 $data_so/SMPD_non_eighted_heidke_skill_individual_realization.py

#cat $data_s/SMPD_heidke_skill.py | sed 's|reanalysis_select|MERRA|g' > $data_so/SMPD_heidke_skill.py && python3 $data_so/SMPD_heidke_skill.py

###################### CASE STUDIES ######################################
reanalysis_list=("NLDAS" "GEFSv12" "MERRA")

for reanalysis in "${reanalysis_list[@]}";do

####################################################  start_date    end_date    weighted RZSM==True   reanalysis

python3 $data_s/CASE_STUDY_with_Albertprojection.py '2012-04-04' '2012-06-27' "True" $reanalysis
python3 $data_s/CASE_STUDY_with_Albertprojection.py '2012-04-04' '2012-06-27' "False" $reanalysis

python3 $data_s/CASE_STUDY_with_Albertprojection.py '2017-04-19' '2017-07-08' "True" $reanalysis
python3 $data_s/CASE_STUDY_with_Albertprojection.py '2017-04-19' '2017-07-08' "False" $reanalysis

python3 $data_s/CASE_STUDY_with_Albertprojection.py '2019-07-15' '2019-10-01' "False" $reanalysis
python3 $data_s/CASE_STUDY_with_Albertprojection.py '2019-07-15' '2019-10-01' "True" $reanalysis;
done



#Looking at RZSM anomalies during CASE STUDY PERIODS

python3 $data_s/RZSM_obs_anomalies_2012_MEM.py
python3 $data_s/RZSM_obs_anomalies_2012_individual_realizations.py

#assessing random probabilistic metrics
python3 $data_s/RZSM_probability_distribution_weighted.py
python3 $data_s/RZSM_probability_distribution_non_weighted.py

#CSI, TPR, and FPR for individual grid cells
python3 $data_s/csi_tpr_fpr_grid_cell.py





################ NOT DOING ANYMORE ##########################################
#ETo all model
#python $data_s/ETo_reference_ET_Penman_ALL_MODEL.py 
#ETo all model
#python $data_s/ETo_reference_ET_Penman_ALL_MODEL_no_NRL.py 

#var_list=("ETo_Penman" "windspeed" )
#This will combine all the models together so that we can make the mean of all models. This may be screwing things up with ETo specifically.
#for var in "${var_list[@]}";do
#cat $data_s/create_multi_model_mean_files.py | sed 's|var_name|'${var}'|g' > $data_so/"$var"_ALL_MODEL_average.py && python3 $data_so/"$var"_ALL_MODEL_average.py;done

#var_list=( "windspeed" "actual_vapor_pressure" "srad")
#for var in "${var_list[@]}";do
#cat $data_s/create_multi_model_mean_files_no_NRL.py | sed 's|var_name|'${var}'|g' > $data_so/"$var"_ALL_MODEL_average_no_NRL.py && python3 $data_so/"$var"_ALL_MODEL_average_no_NRL.py;done



#Now make the multi-model mean for the anomalies, then calculate the anomaly correlation coefficient
#{
#var_list_ETo=("ETo_Penman" "tas" "actual_vapor_pressure" "srad" "pr" "windspeed" "zg")
#for var in "${var_list_ETo[@]}";do
#cat $data_s/anomaly_mean_all_variables_ALL_MODEL.py | sed 's|var_name|'${var}'|g' > $data_so/anomaly_mean_"$var"_ALL_MODEL.py && python3 $data_so/anomaly_mean_"$var"_ALL_MODEL.py
#cat $data_s/MERRA0_reformat_to_SubX_ALL_MODEL.py | sed 's|var_name|'${var}'|g' > $data_so/MERRA0_reformat_to_SubX_"$var"_ALL_MODEL.py && python3 $data_so/MERRA0_reformat_to_SubX_"$var"_ALL_MODEL.py
#cat $data_s/anomaly_correlation_all_variables_ALL_MODEL.py | sed 's|var_name|'${var}'|g' > $data_so/anomaly_correlation_"$var"_ALL_MODEL.py && python3 $data_so/anomaly_correlation_"$var"_ALL_MODEL.py;done
}
################################################END#############################

