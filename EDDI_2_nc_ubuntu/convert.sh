#!/bin/bash

module load cdo

outDir=""

#convert EDDI 1wk .asc to a .nc4 file
for f in *.asc;
do
    year=`echo $f | cut -c16-19`
    month=`echo $f | cut -c20-21`
    day=`echo $f | cut -c22-23`
    outName=`echo $f | cut -c1-23`

    cdo -f nc -settaxis,${year}-${month}-${day},00:00:00,1day -setname,EDDI -input,SMERGE_4_EDDI.grd ${outDir}/${outName}.nc4 < $f;
done
