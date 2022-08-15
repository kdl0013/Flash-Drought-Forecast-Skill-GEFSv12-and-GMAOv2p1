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
processors=8

data_d=$main_directory/Data
data_s=$main_directory/Scripts
mkdir -p $data_d/SubX/GMAO
mkdir -p $data_s 
mkdir -p $data_d/ERA5_land/raw_data/dayTotal
mkdir -p $data_d/EDDI/convert_2_nc

mask=$data_s/CONUS_mask

###############################################################################
#ERA5-Land Soil Moisture and Evaporative Demand Data download, compress
rzsm_obs=$data_d/ERA5_land/raw_data
cd $rzsm_obs

python3 $data_s/00_download_ERA5_Land.py &
sleep 600

#because the file sizes are so large, compress along the way
var_avg=(volumetric*.nc surface*.nc 2m*.nc)
for f in "${var_avg[@]}";
do
if [ ! -f dayTotal/$f ];
then
#sleep 1.5
#Remap to SubX coordinates, take day average, minimze bbox 
cdo -sellonlatbox,235,293,24,50 -remapcon,$mask/CONUS_mask.grd -dayavg $f 1_$f
#compress
ncks -4 -L 1 1_$f compress_$f
rm 1_$f && cp compress_$f dayTotal/$f && mv compress_$f $f
fi;
done
#Get the daily total precipitation, not daily average
for f in total_precipitation*.nc;
do
if [ ! -f dayTotal/$f ];
then
#sleep 1.5
#Remap to SubX coordinates, take day sum (total), minimze bbox 
cdo -sellonlatbox,235,293,24,50 -remapcon,$mask/CONUS_mask.grd -daysum $f 1_$f
#compress
ncks -4 -L 1 1_$f compress_$f
rm 1_$f && cp compress_$f dayTotal/$f && mv compress_$f $f
fi;
done 


#Next, mergetime files to be 1 single file
ulimit -n 8000
cd $rzsm_obs/dayTotal
cdo mergetime surface_net* $data_d/ERA5_land/net_radiation.nc
ncrcat total_pre* $data_d/ERA5_land/precipitation.nc
ncrcat 2m* $data_d/ERA5_land/tavg.nc
cdo mergetime volumetric*layer_1* $data_d/ERA5_land/soil1.nc
cdo mergetime volumetric*layer_2* $data_d/ERA5_land/soil2.nc
cdo mergetime volumetric*layer_3* $data_d/ERA5_land/soil3.nc
cdo mergetime volumetric*layer_4* $data_d/ERA5_land/soil4.nc

###############################################################################
#SubX GEOS-5, FIMr1p1, CCSM4 hindcasts/forecasts, and ESRL forecasts
#return from Cheyenne UCAR cluster
return_SubX_files () {
rsync -Pa klesinger@cheyenne.ucar.edu:/glade/u/home/klesinger/SubX/S* /home/kdl/Insync/OneDrive/NRT_CPC_Internship/Scripts/NCAR_scripts &&
rsync -Pa klesinger@cheyenne.ucar.edu:/glade/u/home/klesinger/SubX/*/compress_resize/* /home/kdl/Insync/OneDrive/NRT_CPC_Internship/Data/SubX/from_easley/
}
return_SubX_files

move_subX=$data_d/SubX
cd $move_subX
models=(GMAO ESRL EMC RSMAS)

for model in "${models[@]}";
do mkdir $model
mv *$model* $model/;
done

ln -s $move_subX/EMC/* $move_subX/GEFSv12 #only because the forecast script was created first

#Priestley-Taylor coefficient
cdo -f nc -sellonlatbox,235,293,24.0,50.0 -remapcon $mask/CONUS_mask.grd -setname,PT_coef -input,$data_s/SMERGE_4_EDDI.grd out.nc4
###############################################################################
#Convert GMAO profile depth to a netcdf file
cd $data_d/SubX/GMAO_soil_moisture_depth
cat prof_depth_1x1.dat | sed -e 's|NaN|0|g' > prof_depth_changed.dat

{

############## Elevation dataset preprocess ###########
mkdir $data_d/elevation && cd $data_d/elevation
data_e=$data_d/elevation

cdo -sellonlatbox,235,293,24,50 -remapcon,$data_s/CONUS_mask.grd $data_e/elev.1-deg.nc $data_e/elev_regrid.nc

############## Evaporative Demand Drought Index ##########################
{
echo "Working on EDDI pre-processing"

eddi=$data_d/EDDI
mkdir ${eddi}/convert_2_nc
cd ${eddi}

cat $data_s/0c_download_EDDI.py | sed 's|main_dir|'${main_directory}'|g' \
> $data_s/TMP_0c_download_EDDI.py 

#Create wget script for EDDI files
python3 $data_s/TMP_0c_download_EDDI.py
rm $data_s/TMP_0c_download_EDDI.py

bash $data_s/0d_download_EDDI_wget.sh
}
#Some downloads produce double on my computer, this moves the doubles so you
#could verify if the actual file is present

#for f in *.asc; 
#do  
#    if [[ ${#f} -gt 28 ]]; 
#then 
#    rm "$f"; fi;  done
 
#convert to a .nc4 file
for f in *.asc;
do
    year=`echo $f | cut -c16-19`
    month=`echo $f | cut -c20-21`
    day=`echo $f | cut -c22-23`
    outName=`echo $f | cut -c1-23`

       cdo -f nc -sellonlatbox,235,293,24.0,50.0 -remapcon,$data_s/CONUS_mask.grd -settaxis,${year}-${month}-${day},00:00:00,1day -setname,EDDI -input,$data_s/SMERGE_4_EDDI.grd convert_2_nc/${outName}.nc4 < $f;
done

ulimit -n 7000
#Merge files

cdo mergetime convert_2_nc/*.nc4 convert_2_nc/eddi_merged.nc4

for f in *.asc;
do
    year=`echo $f | cut -c1-4`
    month=`echo $f | cut -c5-6`
    day=`echo $f | cut -c7-8`
    outName=`echo $f | cut -c1-8`

       cdo -f nc -remapcon,$data_s/CONUS_mask.grd -settaxis,${year}-${month}-${day},00:00:00,1day -setname,EDDI -input,$data_s/SMERGE_4_EDDI.grd EDDI_nc_file/${outName}.nc4 < $f;
done

#Download GMAO GEOS-5 year 2017-2021
python3 $data_s/00_create_wget_GEOS_2017_2021_years.py
#save a wget_GMAO_2021.sh script. I moved this into the download location because 
#wget parameter (-O) wasn't saving properly
cd /home/kdl/Insync/OneDrive/NRT_CPC_Internship/Data/SubX/GMAO/preprocess_2017_2021
bash wget_GMAO_2021.sh

#Now shrink files by converting to smaller grid.
for f in *.nc;do
cdo -sellonlatbox,235,293,24,50 -remapcon,$mask/CONUS_mask.grd $f "$f"4 &&
mv "$f"4 $f;
done

#OLD CODE, SMERGE only goes through year 2019, switching to GLEAM soil moisture instead
################### SMERGE soil moisture ##########################

#cd $data_s
#cat 0a_download_SMERGE_RZSM.py | sed 's|main_dir|'${main_directory}'|g' | \
#sed 's|procs|'${processors}'|g' > TMP_Oa_download_SMERGE_RZSM.py

#echo Starting SMERGE Soil Moisture Download
#Run file

#python3 TMP_Oa_download_SMERGE_RZSM.py
  









#OLD SMERGE RZSM
{
cd $data_d/SMERGE_SM/Raw_data
#Remove files that had wildcard expressions that were doubly downloaded

ulimit -n 9500
cdo mergetime *.nc4 smerge_sm_merged.nc4
cdo -sellonlatbox,235,293,24,50 -remapcon,$data_s/CONUS_mask.grd smerge_sm_merged.nc4 smerge_sm_merged_remap.nc4
}



############## ONLY FOR OLD GLEAM CODE ######################33

#Shrink files first, merging them makes a very large dataset
for f in *.nc;do
cdo -sellonlatbox,235,293,24,50 -remapcon,$mask/CONUS_mask_GLEAM.grd $f a_$f;
done


cdo mergetime a_SMroot* SMroot_merged.nc4
cdo mergetime a_SMsurf* SMsurf_merged.nc4

#TODO: Only run this command once. The next time, it will append more data that makes 
#the file unreadable. If so, just change the lon,a,c to lon,o,c for all variables (and for lat)
for f in *.nc;do
ncatted -a standard_name,lon,a,c,"longitude" $f
ncatted -a long_name,lon,a,c,"longitude" $f
ncatted -a units,lon,a,c,"degrees_east" $f
ncatted -a axis,lon,a,c,"lon" $f

ncatted -a standard_name,lat,a,c,"latitude" $f
ncatted -a long_name,lat,a,c,"latitude" $f
ncatted -a units,lat,a,c,"degrees_north" $f
ncatted -a axis,lat,a,c,"lat" $f;
done



#TEst convert preiestley
 cdo -f nc -sellonlatbox,235,293,24.0,50.0 -settaxis,2000-01-01,00:00:00,1day -setname,EDDI -input,$mask/Priestley.grd test.nc4 < pt_coeff_ascii.txt
 
nc

#Move into the directory where files will be saved
cd $data_d/gridMET

#Download METDATA files
bash $data_s/0b_download_metdata_wget.sh
#Download year 2021 seperately from website https://www.northwestknowledge.net/metdata/data/

#Merge together all years of each variable, remap areal conservative, select CONUS
for var in srad tmmn tmmx vs rmin rmax;
do
echo Starting variable $var
 cdo -sellonlatbox,235,293,24,50 -remapcon,$mask/CONUS_mask.grd -mergetime ${var}_*.nc ${var}_remap_final.nc;
done

cdo ensavg tmmn_remap_final.nc tmmx_remap_final.nc tavg_remap_final.nc

echo gridMET files are downloaded and pre-processed
}
