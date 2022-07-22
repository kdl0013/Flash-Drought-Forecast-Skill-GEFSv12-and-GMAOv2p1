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
model=GMAO 
main_directory='/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
processors=7

data_d=$main_directory/Data
data_s=$main_directory/Scripts

subx=$data_d/SubX/${model}
mkdir $subx/anomaly
#mkdir $data_s/cython_scripts


# --- Functions
make_reference_ET () {
name=${1::-3}
cat $data_s/$1 | sed 's|main_dir|'${main_directory}'|g' | 
    sed 's|procs|'${processors}'|g' | sed 's|model_name|'${model}'|g' > $data_s/TMP_${name}_${model}.py
    
python3 $data_s/TMP_${name}_${model}.py 
}

make_reference_ET "1b1_make_reference_ET.py"



#convert soil moisture to m3/m3
convert_RZSM () {
cat 1b2_convert_RZSM_to_m3_m3.py | sed 's|main_dir|'${main_directory}'|g' > $data_s/TMP_1b2_convert_RZSM_to_m3_m3.py

python3 $data_s/TMP_1b2_convert_RZSM_to_m3_m3.py
}

convert_RZSM

#Remove uneeded S dimension
remove_S_dimension () {
cat 1b3_remove_S_dimsension.py | sed 's|model_name|'${model}'|g' | sed 's|main_dir|'${main_directory}'|g'> $data_s/TMP_1b3_remove_S_dimsension.py 

python3 TMP_1b3_remove_S_dimsension.py
}

remove_S_dimension

make_empty_anomaly_files () {
name=${1::-3}
var=$2

cat $data_s/$1 | sed 's|main_dir|'${main_directory}'|g' | 
    sed 's|procs|'${processors}'|g' | sed 's|model_name|'${model}'|g' | sed 's|variables|'${var}'|g' > $data_s/TMP_${name}_${model}.py
    
python3 $data_s/TMP_${name}_${model}.py 
}

make_empty_anomaly_files "1b4_make_empty_anomaly_files.py" "ETo"
make_empty_anomaly_files "1b4_make_empty_anomaly_files.py" "RZSM"
make_empty_anomaly_files "1b4_make_empty_anomaly_files.py" "EDDI"

######## Reference ETo, EDDI, gridMET data creation. Historical and SubX #######
#Create multiple scripts for multiprocessing of EDDI, RZSM, and ETo anomalies
#Can also insert the total number of models. $1 = name of file, $2 = 3 models for GMAO

#Changes the initialization start date for processing several files
#Multiprocessing in python does not work with this data (CPU is constrained to 12 files)
###I can do at most 3 files at a time due to how the files need to be accessed up to 3 months lead/lag
###Using 3 files at a time makes it too slow to compute on my personal computer, I could split it up

mk_anomaly_script_RZSM_ETo_EDDI () {
date_arr=(0 25 50)
step=25
name=${1::-3}

for val in "${date_arr[@]}";do
    for mod in {0..3};do
cat $1 | sed 's|main_dir|'${main_directory}'|g' | sed 's|start_init|'${val}'|g' | sed 's|init_step|'${step}'|g' | sed 's|model_name|'${model}'|g' | sed 's|variable_|'$2'|g' > TMP_${name}_"$2"_step"$val"_${model}.py;
    done;
done
}

#Create tmp scripts to run later
#mk_anomaly_script_RZSM_ETo_EDDI '1c_EDDI.pyx'
mk_anomaly_script_RZSM_ETo_EDDI "1c_make_anomaly.py" "RZSM"
mk_anomaly_script_RZSM_ETo_EDDI "1c_make_anomaly.py" "ETo"


#Run all three steps at a single time. Takes 30GB of RAM to run 3 at a time.
run_anomaly_RZSM_ETo_EDDI () {
cd $data_s
for file in TMP_1c_make_anomaly_$1*.py;do
python3 $file &
done
}

run_anomaly_RZSM_ETo_EDDI "RZSM"
run_anomaly_RZSM_ETo_EDDI "ETo"


#Insert EDDI here

#Check for missing anomaly data in files
check_missing_anomaly_data () {
cat 1e_check_for_missing_data.py | sed 's|main_dir|'${main_directory}'|g' | sed 's|model_name|'${model}'|g' > TMP_1e_check_for_missing_data.py

python3 TMP_1e_check_for_missing_data.py
}
check_missing_anomaly_data

#Make gridMET and ETo anomalies using same methodology as subx
make_obs_anomaly () {
cat 1f_make_observation_anomalies.py | sed 's|main_dir|'${main_directory}'|g' > TMP_1f_make_observation_anomalies.py

python3 TMP_1f_make_observation_anomalies.py
}
make_obs_anomaly



#Create new dataset to compare the skill between models
reformat_observations_to_SubX_format() {
cat 1i_Obs_reformatted_to_SubX.py | sed 's|main_dir|'${main_directory}'|g' | sed 's|procs|'${processors}'|g' > TMP_1i_Obs_reformatted_to_SubX.py

python3 TMP_1i_Obs_reformatted_to_SubX.py
}

reformat_observations_to_SubX_format #call function

#convert percentiles to SMPD binary occurence of flash drought
make_SMPD_FD_classification () {
name=${1::-3}
cat $data_s/$1 | sed 's|main_dir|'${main_directory}'|g' | 
    sed 's|procs|'${processors}'|g' | sed 's|model_name|'${model}'|g' > $data_s/TMP_${name}_${model}.py
    
python3 $data_s/TMP_${name}_${model}.py 
}

make_SMPD_FD_classification "1n_RZSM_FD_classification_SubX.py"













make_mean_for_acc () {
date_arr=(0 25 50)
step=25
name=${1::-3}

for val in "${date_arr[@]}";do
    for mod in {0..3};do
cat $1 | sed 's|main_dir|'${main_directory}'|g' | sed 's|start_init|'${val}'|g' | sed 's|init_step|'${step}'|g' | sed 's|model_name|'${model}'|g' > TMP_${name}_step"$val"_${model}.py;
    done;
done

}

make_mean_for_acc "1c2_anomaly_RZSM_mean_for_ACC.py"

#chmod +x *.pyx
chmod +x *.py

run_anomaly_RZSM_ETo_EDDI () {
cd $data_s
for file in TMP_*anomaly_$1*.py;do
python3 $file &
done
}

run_anomaly_RZSM_ETo_EDDI "RZSM"

run_mean_RZSM_ETo_EDDI () {
cd $data_s
for file in TMP_*anomaly_$1*mean*.py;do
python3 $file &
done
}

run_mean_RZSM_ETo_EDDI "RZSM"
#Call cython_setup.py, insert file name, and compile 
#Insert only name of variable
compile_cython () {
for file in TMP*"$1"*.pyx;
do
cat $data_s/cython_setup.py | sed 's|insert_cython_file|'$file'|g' > cy_${file}
python3 cy_${file} build_ext --inplace
rm cy_${file};
done
#copy all created .so and .c files into a new folder to keep directory smaller
mv *.so cython_scripts/
mv *.c cython_scripts/
}

#Compile code using cython
compile_cython 'ETo' &
#compile_cython 'EDDI' &
compile_cython 'RZSM'

#Create a new python script that imports the newly created .c extension
#The file starts running once the module is imported
create_and_run_cython_scripts () {
cat << EOF > TMP_run_$1.py
#/usr/bin/env python3

import $1 
exit()
EOF
python3 TMP_run_$1.py
}


run_anomaly_RZSM_ETo_EDDI () {
cd $data_s/cython_scripts
for file in TMP_*$1*.c;do
create_and_run_cython_scripts $file &
done

cd $data_s
}


run_anomaly_RZSM_ETo_EDDI "ETo" 
run_anomaly_RZSM_ETo_EDDI "RZSM"
#run_anomaly_RZSM_ETo_EDDI "EDDI"




#Fix any broken SubX anomaly files
fix_broken_subX_anomalies () {
cat $1 | sed 's|main_dir|'${main_directory}'|g' | sed 's|model_name|'${model}'|g' > TMP_"$model"_$1

python3 TMP_"$model"_$1
}

fix_broken_subX_anomalies "1f_fix_individual_anomalies_ETo.py"
fix_broken_subX_anomalies "1g_fix_individual_anomalies_RZSM.py"


#Create gridMET and SMERGE anomlies (observations)
make_obs_anomalies() {
cat 1h_make_gridMET_anomalies.py | sed 's|main_dir|'${main_directory}'|g' > TMP_1h_make_gridMET_anomalies.py

cat 1h2_make_gridMET_anomalies_fix_january.py | sed 's|main_dir|'${main_directory}'|g' > TMP_1h2_make_gridMET_anomalies_fix_january.py

python3 TMP_1h_make_gridMET_anomalies.py
python3 TMP_1h2_make_gridMET_anomalies_fix_january.py
}

make_obs_anomalies

















############## Lagged Average Ensemble ##############
cat $data_s/3a_lagged_average_ensemble.py | sed 's|main_dir|'${main_directory}'|g' > $data_s/3a_TMP_lagged_average_ensemble.py

python3 $data_s/3a_TMP_lagged_average_ensemble.py
rm $data_s/3a_TMP_lagged_average_ensemble.py



#EVERYTHING BELOW THIS LINE IS EXTRA CODE

#Run cython file for each .c file (this doesnt' allow loop to continue)
run_cythonB () {
for file in TMP_*$1*.c;
do
run_cythonA $file &
done
}

#Test 2 (probably won't work)
run_anomaly () {
x=1
while [ x -le 24 ];do
for file in TMP_1e_anomaly_RZSM_step*_${model}.pyx;do
python3 $file &
wait; done
x=$(( x++ ));
done
}





#Copy already completed files back into correct directory, this works if the 
#anomaly calculation was halted (which would break some files while saving)
copy_broken_files () {
for mod in {0..3};do
cp ${subx}/EDDI_mod"$mod"/already_completed/*.npy ${subx}/EDDI_mod"$mod"
cp ${subx}/RZSM_anomaly_mod"$mod"/already_completed/*.npy ${subx}/RZSM_anomaly_mod"$mod"
rm ${subx}/ETo_anomaly_mod"$mod"/*.npy
cp ${subx}/ETo_anomaly_mod"$mod"/already_completed/*.npy ${subx}/ETo_anomaly_mod"$mod";
done
}

copy_broken_files


#Copy empty .npy files created in 1a_ETo.py into appropriate folders
cp_npy_files () {
for mod in {0..3};do
mkdir -p ${subx}/$1_mod${mod}/already_completed
cp ${subx}/$1*.npy ${subx}/$1_mod${mod};
done
}


#Copy empty .npy files before processing
cp_npy_files "ETo_anomaly"
cp_npy_files "EDDI"
cp_npy_files "RZSM_anomaly"





EDDI_function () {
for file in TMP_1b_EDDI_step*_${model}.py;do 
python3 $file &
done
}

for i in {1..24};do
EDDI_function
wait;
done


#If EDDI_function breaks when no fully completed, run copy_broken_files and restart EDDI_function
copy_broken_files () {
for mod in {0..3};do
cp ${subx}/EDDI_mod"$mod"/already_completed/*.npy ${subx}/EDDI_mod"$mod";
done
}


#Combine all EDDI model outputs
cat 1c_EDDI_stitch_models.py | sed 's|main_dir|'${main_directory}'|g' | sed 's|mod_name|'${model}'|g' > TMP_1c_EDDI_stitch_models_${model}.py 

python3 TMP_1c_EDDI_stitch_models_${model}.py

############## Reference ETo,EDDI,RZSM scatterplot  ###########
cat $data_s/1d_RZSM_and_scatterplots.py | sed 's|main_dir|'${main_directory}'|g' | \
    sed 's|procs|'${processors}'|g' | sed 's|mod_name|'${model}'|g' > $data_s/TMP_1d_RZSM_and_scatterplots_${model}.py 

python3 $data_s/TMP_1d_RZSM_and_scatterplots_${model}.py


#####Next, process anomaly EDDI and RZSM, make each model into it's own script for speed
#Dont do anomaly for EDDI. It is already standardized.
#for i in {0..3};do
#mkdir ${subx}/EDDI_anomaly_mod${i}
#cp ${subx}/EDDI_anomaly_*.npy ${subx}/EDDI_anomaly_mod${i};
#done
#rm ${subx}/EDDI_anomaly_*.npy

#Make new directories and copy files to begin multiprocessing on RZSM anomaly
for i in {0..3};do
mkdir -p ${subx}/RZSM_anomaly_mod${i}/already_completed
cp ${subx}/RZSM_anomaly_*.npy ${subx}/RZSM_anomaly_mod${i};
done
#rm ${subx}/RZSM_anomaly_*.npy

for val in "${date_arr[@]}";do
    for mod in {0..3};do
cat 1e_anomaly_RZSM.py | sed 's|main_dir|'${main_directory}'|g' | sed 's|start_init|'${val}'|g' | sed 's|init_step|'${step}'|g' | sed 's|model_number|'${mod}'|g' | sed 's|model_name|'${model}'|g' > TMP_1e_anomaly_RZSM_step"$val"_mod"$mod"_${model}.py;
    done;
done

{


for i in {1..24};do
RZSM_anomaly_function
wait;
done
}


x=1
while [ x -le 24 ];
do
for file in TMP_1e_anomaly_RZSM_step*_${model}.py;
do
python3 $file &
wait; done
x=$(( x++ ));
done

copy_broken_RZSM () {
for mod in {0..3};do
cp ${subx}/RZSM_anomaly_mod"$mod"/already_completed/*.npy ${subx}/RZSM_anomaly_mod"$mod";
done
}

#TODO: Import back the data from HPC cluster and stitch together
#current issue, HPC is very slow...not sure why

for i in {0..3};do
mkdir -p ${subx}/ETo_anomaly_mod${i}/already_completed
cp ${subx}/ETo_anomaly_*.npy ${subx}/ETo_anomaly_mod${i};
done

#HPC date steps allowed
date_arr=(0 25 50)
step=25
for val in "${date_arr[@]}";do
    for mod in {0..3};do
cat 1f_anomaly_ETo.pyx | sed 's|main_dir|'${main_directory}'|g' | sed 's|start_init|'${val}'|g' | sed 's|init_step|'${step}'|g' | sed 's|model_number|'${mod}'|g' | sed 's|model_name|'${model}'|g' > TMP_1e_anomaly_ETo_step"$val"_mod"$mod"_${model}.pyx;
    done;
done

ETo_anomaly_function () {
for file in TMP_1f_anomaly_ETo_step*_${model}.py;do 
python3 $file &
done
}

for i in {1..24};do
ETo_anomaly_function
wait;
done
}

copy_broken_ETo () {
for mod in {0..3};do
cp ${subx}/ETo_anomaly_mod"$mod"/already_completed/*.npy ${subx}/ETo_anomaly_mod"$mod";
done
}
####Be careful running code below
#rm ${subx}/EDDI* -r
#rm ${data_s}/TMP_*
import os
from glob import glob
import numpy as np

for f in sorted(glob('ETo*.npy')):
    print(np.load(f,allow_pickle=True))

ff = []
for f in sorted(glob('RZSM*.npy')):
    try:
        np.load(f,allow_pickle=True)
    except ValueError:
        ff.append(f)



#Can delete this one
run_anomaly () {
cd $data_s/cython_scripts

x=1
while [ $x -le 24 ];do
for file in TMP_*$1*.c;do
create_and_run_cython_scripts $file &
done; wait
x=$(( x++ ));
done
}

