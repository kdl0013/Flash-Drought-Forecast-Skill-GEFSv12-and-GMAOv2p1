#!/bin/bash

#Download SMERGE soil moisture, gridMET PET, EDDI evaporative demand, and elevation
#To download SubX data, use script download_SubX_4_models.sh (run in HPC)

#My own environment is spyder
conda activate spyder 
module load cdo
module load nco

#MUST SET inputs: main_directory and number of processors for CPU parallelism

#You can make main_directory path any path you would like and this will get
#all other directories in order. 

main_directory='/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
processors=8

data_d=$main_directory/Data
data_s=$main_directory/Scripts
mkdir -p $data_d/SubX/GMAO
mkdir -p $data_s 
mkdir -p $data_d/SMERGE_SM/Raw_data
mkdir -p $data_d/gridMET/ETo_SubX_values
mkdir -p $data_d/EDDI/convert_2_nc
################### SMERGE soil moisture ##########################
{
cd $data_s
cat 0a_download_SMERGE_RZSM.py | sed 's|main_dir|'${main_directory}'|g' | \
sed 's|procs|'${processors}'|g' > 0a_TMP_download_SMERGE_RZSM.py

echo Starting SMERGE Soil Moisture Download
#Run file

python3 0a_TMP_download_SMERGE_RZSM.py
rm 0a_TMP_download_SMERGE_RZSM.py
}  
#####Remap RZSM files to a 1 x 1 grid and then select only CONUS
#Permutate dimensions to be the same as SubX
{
cd $data_d/SMERGE_SM/Raw_data
#Remove files that had wildcard expressions that were doubly downloaded
mkdir wild_card_files

for f in *.nc4; 
do  
    if [[ ${#f} -gt 12 ]]; 
then 
    mv "$f" wild_card_files/; fi;  done

ulimit -n 6500
cdo mergetime *.nc4 smerge_sm_merged.nc4
cdo -sellonlatbox,235,293,24,50 -remapcon,$data_s/regrid_CONUSmask_all_variables.grd smerge_sm_merged.nc4 smerge_sm_merged_remap.nc4
}
################### gridMET METDATA ##########################
echo "Starting METDATA download of PET (alfalfa)"
#Move into the directory where files will be saved
cd $data_d/gridMET
mkdir remap_merged

#Download METDATA files
bash $data_s/0b_download_metdata_wget.sh

#Merge together all years of each variable, remap areal conservative, select CONUS
for var in srad tmmn tmmx vs rmin rmax;
do
echo Starting variable $var
 cdo -sellonlatbox,235,293,24,50 -remapcon,$data_s/regrid_CONUSmask_all_variables.grd -mergetime ${var}_*.nc ${var}_remap_final.nc;
done

cdo ensavg tmmn_remap_final.nc tmmx_remap_final.nc tavg_remap_final.nc

echo gridMET files are downloaded and pre-processed
############## Elevation dataset preprocess ###########
mkdir $data_d/elevation && cd $data_d/elevation
data_e=$data_d/elevation

cdo -sellonlatbox,235,293,24.0,50.0 -remapcon,$data_s/regrid_CONUSmask_all_variables.grd $data_e/elev.1-deg.nc $data_e/elev_regrid.nc

############## Evaporative Demand Drought Index ##########################
echo "Working on EDDI pre-processing"

eddi=$data_d/EDDI
cd eddi

cat $data_s/0c_download_EDDI.py | sed 's|main_dir|'${main_directory}'|g' \
> $data_s/0c_TMP_download_EDDI.py 

#Download EDDI files
python3 $data_s/0c_TMP_download_EDDI.py
rm $data_s/0c_TMP_download_EDDI.py

bash $data_s/0d_download_EDDI_wget.sh

mkdir wild_card_files

#Some downloads produce double on my computer, this moves the doubles so you
#could verify if the actual file is present
for f in *.asc; 
do  
    if [[ ${#f} -gt 27 ]]; 
then 
    mv "$f" wild_card_files/; fi;  done

#convert to a .nc4 file
for f in *.asc;
do
    year=`echo $f | cut -c16-19`
    month=`echo $f | cut -c20-21`
    day=`echo $f | cut -c22-23`
    outName=`echo $f | cut -c1-23`

       cdo -f nc -sellonlatbox,235,293,24.0,50.0 -remapcon,$data_s/regrid_CONUSmask_all_variables.grd -settaxis,${year}-${month}-${day},00:00:00,1day -setname,EDDI -input,$data_s/SMERGE_4_EDDI.grd convert_2_nc/${outName}.nc4 < $f;
done

#Merge files
cd ${eddi}/convert_2_nc
cdo mergetime *.nc4 eddi_merged.nc4































