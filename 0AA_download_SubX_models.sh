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
model_array=(GMAO ESRL RSMAS EMC)


#Move data to HPC (easley)
rsync -Pa  ~/Insync/OneDrive/NRT_CPC_Internship/ kdl0013@easley.auburn.edu:/home/kdl0013/NRT_CPC_Internship


#Create wget scripts for each variable
#TODO:Return netcdf file forecasts from GEFSv12 from EASLEY####################
echo Mi
return_easley_files () {
#Return scripts for downloading and converting data
rsync -Pa kdl0013@easley.auburn.edu:/home/kdl0013/GEFSv12/* ~/Insync/OneDrive/NRT_CPC_Internship/Data/SubX/*/raw_ensemble_files 
rsync -Pa kdl0013@easley.auburn.edu:/home/kdl0013/process_GEFS/* ~/Insync/OneDrive/NRT_CPC_Internship/Scripts/EASLEY_HPC/ && \
rsync -Pa kdl0013@easley.auburn.edu:/home/kdl0013/wget* ~/Insync/OneDrive/NRT_CPC_Internship/Scripts/EASLEY_HPC/
}
#Contains years 2020 - present from easley
return_easley_files #for GEFSv12 only

#maybe remove d35_dswrf_sfc_2018072500_p04.nc4          tmp_2m/d35_tmp_2m_2003031900_p08.nc4

echo An
#TODO:Return netcdf file forecasts from CASPER (all models, part of EMC)#######
return_CASPER_SUBx_files () {
mkdir -r $1
mkdir -r $2
rsync -Pa klesinger@cheyenne.ucar.edu:/glade/u/home/klesinger/SubX/ECCC/compress_resize/* $2
rsync -Pa klesinger@cheyenne.ucar.edu:/glade/u/home/klesinger/SubX/S* $1
rsync -Pa klesinger@cheyenne.ucar.edu:/glade/u/home/klesinger/SubX/wget* $1
}
return_CASPER_SUBx_files $scripts $outData

rsync -Pa klesinger@cheyenne.ucar.edu:/glade/u/home/klesinger/SubX/ECCC/compress_resize/* /home/kdl/Insync/OneDrive/NRT_CPC_Internship/Data/SubX/fromCasper/ECCC




remove_s_dim () {
#TODO: Remove s dimension for later processing
rm -r $outData/S_dim_removed/*$1*
#ncks operators don't like overwriting in python without needing input from user
cat 01e_remove_S_dimsension_all_models.py | sed 's|model_name|'${1}'|g'> $data_s/$1_01e_remove_S_dimsension_all_models.py && python3 $data_s/$1_01e_remove_S_dimsension_all_models.py 
rm $outData/S_dim_removed/a_* $outData/S_dim_removed/b_*

}

remove_s_dim "ECCC"




#TODO: Preprocess files sent from ESRL FIMr1p1 agency #########################
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

#rename files
python3 $data_s/ESRL1_rename_files.py
################################################################################


###Preprocess data
#TODO: Make a script to only get the correct days for RSMAS#####################
python3 $data_s/01d_select_RSMAS_dates.py



#TODO: Move files into seperate directories
model_array=(EMC)
for model in "${model_array[@]}";
do mkdir $subx/$model
cp $outData/S_dim_removed/*$model*.nc4 $subx/$model/;
done




python3 $data_s/EMC0_day_average_GEFSv12_HPC.py
python3 $data_s/EMC0_day_sum_GEFSv12_HPC.py
python3 $data_s/EMC1_merge_3_soil_moisture_fields_GEFSv12.py
mv $data_d/SubX/EMC/raw_ensemble_files /media/kdl/Seagate_1/CPC_project/raw_ensemble_files #to save on storage space
mv $data_d/SubX/EMC/raw_ensemble_files /media/kdl/Seagate_1/CPC_project/raw_ensemble_files2/*


make_ETo_Penman_single () {
cat ETo_reference_ET_Penman.py | sed 's|model_name|'${1}'|g'> $data_s/"$1"_ETo_reference_ET_Penman.py && python3 $data_s/"$1"_ETo_reference_ET_Penman.py
}
make_ETo_Penman_single "GMAO"


make_ETo_Penman () {
model_array=$1
for model in "${model_array[@]}";
do cat ETo_reference_ET_Penman.py | sed 's|model_name|'${model}'|g'> $data_s/"$model"_ETo_reference_ET_Penman.py && python3 $data_s/"$model"_ETo_reference_ET_Penman.py;
done
}
make_ETo_Penman $model_array


#Make ETo mean 
ETo_mean_anomaly_and_correlation () {
#create the mean file (to subtract and make anomaly)
cat $data_s/anomaly_mean_ETo_not_EMC.py | sed 's|model_name|'${1}'|g' | sed 's|evap_equation|'${2}'|g' > $data_s/$1_anomaly_mean_$2_not_EMC.py && python3 $data_s/$1_anomaly_mean_$2_not_EMC.py
#Make the observations in same format to go with anomaly correlation
#cat $data_s/MERRA0_reformat_to_SubX_ETo.py | sed 's|model_name|'${1}'|g' | sed 's|evap_equation|'${2}'|g' > $data_s/$1_MERRA0_reformat_to_SubX_ETo_$2.py && python3 $data_s/$1_MERRA0_reformat_to_SubX_ETo_$2.py
#Create anomaly 
#cat $data_s/anomaly_compute_from_mean_ETo.py | sed 's|model_name|'${1}'|g' | sed 's|evap_equation|'${2}'|g' > $data_s/$1_anomaly_compute_from_mean_ETo_$2.py && python3 $data_s/$1_anomaly_compute_from_mean_ETo_$2.py
#Plot correlation
#cat $data_s/anomaly_correlation_ETo.py | sed 's|model_name|'${1}'|g' | sed 's|evap_equation|'${2}'|g' > $data_s/$1_anomaly_correlation_$2.py && python3 $data_s/$1_anomaly_correlation_$2.py
cat $data_s/CRPS_skill_ETo.py | sed 's|model_name|'${1}'|g' | sed 's|evap_equation|'${2}'|g' > $data_s/$1_CRPS_skill_ETo_$2.py && python3 $data_s/$1_CRPS_skill_ETo_$2.py
#cat $data_s/ACC_skill_ETo.py | sed 's|model_name|'${1}'|g' | sed 's|evap_equation|'${2}'|g' > $data_s/$1_ACC_skill_ETo_$2.py && python3 $data_s/$1_ACC_skill_ETo_$2.py
}

ETo_mean_anomaly_and_correlation "ESRL" "Priestley"
ETo_mean_anomaly_and_correlation "RSMAS" "Priestley"

ETo_mean_anomaly_and_correlation "GMAO" "Priestley"
ETo_mean_anomaly_and_correlation "GMAO" "Penman"







EMC_anomaly_mean () {
#Make ETo mean file
cat $data_s/anomaly_mean_ETo_only_EMC_11_models.py | sed 's|model_name|'${1}'|g' | sed 's|evap_equation|'${2}'|g' > $data_s/$1_anomaly_mean_ETo_$2_11_models.py && python3 $data_s/$1_anomaly_mean_ETo_$2_11_models.py
#Make RZSM mean file
cat $data_s/anomaly_mean_RZSM_only_EMC_11_models.py | sed 's|model_name|'${1}'|g' > $data_s/$1_anomaly_mean_RZSM_11_models.py && python3 $data_s/$1_anomaly_mean_RZSM_11_models.py

#ETo
#cat $data_s/MERRA0_reformat_to_SubX_ETo.py | sed 's|model_name|'${1}'|g' | sed 's|evap_equation|'${2}'|g' > $data_s/$1_MERRA0_reformat_to_SubX_ETo_$2.py && python3 $data_s/$1_MERRA0_reformat_to_SubX_ETo_$2.py
#cat $data_s/anomaly_compute_from_mean_ETo.py | sed 's|model_name|'${1}'|g' | sed 's|evap_equation|'${2}'|g' > $data_s/$1_anomaly_compute_from_mean_ETo_$2.py && python3 $data_s/$1_anomaly_compute_from_mean_ETo_$2.py

#RZSM
#cat $data_s/MERRA0_reformat_to_SubX_RZSM.py | sed 's|model_name|'${1}'|g' > $data_s/$1_MERRA0_reformat_to_SubX_RZSM.py && python3 $data_s/$1_MERRA0_reformat_to_SubX_RZSM.py
#cat $data_s/anomaly_compute_from_mean_RZSM.py | sed 's|model_name|'${1}'|g' > $data_s/$1_anomaly_compute_from_mean_RZSM.py && python3 $data_s/$1_anomaly_compute_from_mean_RZSM.py

#cat $data_s/anomaly_correlation_ETo.py | sed 's|model_name|'${1}'|g' | sed 's|evap_equation|'${2}'|g' > $data_s/$1_anomaly_correlation_$2.py && python3 $data_s/$1_anomaly_correlation_$2.py

cat $data_s/anomaly_correlation_RZSM.py | sed 's|model_name|'${1}'|g'  > $data_s/$1_anomaly_correlation_RZSM.py && python3 $data_s/$1_anomaly_correlation_RZSM.py


}
EMC_anomaly_mean "EMC" "Priestley"

EMC_anomaly_mean () {
#cat $data_s/ACC_skill_ETo.py | sed 's|model_name|'${1}'|g' | sed 's|evap_equation|'${2}'|g' > $data_s/$1_ACC_skill_ETo_$2.py && python3 $data_s/$1_ACC_skill_ETo_$2.py
#cat $data_s/ACC_skill_RZSM.py | sed 's|model_name|'${1}'|g' > $data_s/$1_ACC_skill_RZSM.py && python3 $data_s/$1_ACC_skill_RZSM.py
#cat $data_s/CRPS_skill_ETo.py | sed 's|model_name|'${1}'|g' | sed 's|evap_equation|'${2}'|g' > $data_s/$1_CRPS_skill_ETo_$2.py && python3 $data_s/$1_CRPS_skill_ETo_$2.py
cat $data_s/CRPS_skill_RZSM.py | sed 's|model_name|'${1}'|g' > $data_s/$1_CRPS_skill_RZSM.py && python3 $data_s/$1_CRPS_skill_RZSM.py
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

