#!/bin/sh
# Generate individual files to allow downloading of each ensemble member in parallel.
#
# Created by Ray Bell (https://github.com/raybellwaves).

# User defined variables:
# See http://iridl.ldeo.columbia.edu/SOURCES/.Models/.SubX/ for a list of models availible.
# See http://cola.gmu.edu/kpegion/subx/data/priority1.html or
#http://cola.gmu.edu/kpegion/subx/data/priority2.html for a list of var abbreviations.
# See http://cola.gmu.edu/kpegion/subx/docs/SubXDataQuickReferenceGuide.pdf
#for notes on what presseure level is associated with the data.
#

#Must select the correct conda environment (this is my specific env)
#conda activate spyder

ftype=hindcast # hindcast, forecast
mod=GEFS # 30LCESM1, 46LCESM1, CCSM4, CFSv2, FIMr1p1, GEFS, GEM, GEOS_V2p1, NESM
inst=EMC # CESM,    CESM,     RSMAS, NCEP,  ESRL,    EMC,  ECCC, GMAO,     NRL
var=soilw1 # pr, tas, ts, rlut, ua, va, zg
plev=None # 200, 500, 850, 2m, sfc, toa, None
outdir=/home/kdl/Insync/OneDrive/NRT_CPC_Internship/Data

mkdir $outdir/${inst}_SubX
outdir=$outdir/${inst}_SubX

# Default variables
url=http://iridl.ldeo.columbia.edu/SOURCES/.Models/.SubX/

# Remove any files previously created
rm -rf tmp.py
rm -rf *_e*.py

# Find out the first and last ensemble associated with the model:

# Write python script
cat > tmp.py << EOF
#!/usr/bin/env python
import xarray as xr
_rd = xr.open_dataarray('${url}.${inst}/.${mod}/.${ftype}/.${var}/dods')
print(int(_rd.M.values[0]))
EOF

# Run python script and return the first ensemble
fen=`python tmp.py`
rm -rf tmp.py

# Write python script
cat > tmp.py << EOF
#!/usr/bin/env python
import xarray as xr
_rd = xr.open_dataarray('${url}.${inst}/.${mod}/.${ftype}/.${var}/dods')
print(int(_rd.M.values[-1]))
EOF

# Run python script and return the last ensemble
len=`python tmp.py`
rm -rf tmp.py

for ens in {0..10};do
    # Replace text in python template file for each ensemble member
    cat getSubXdatafull_template.py\
    | sed 's|url|'${url}'|g'\
    | sed 's|outdir|'${outdir}'|g'\
    | sed 's/ftype/'${ftype}'/g'\
    | sed 's/mod/'${mod}'/g'\
    | sed 's/inst/'${inst}'/g'\
    | sed 's/var/'${var}'/g'\
    | sed 's/plev/'${plev}'/g'\
    | sed 's/ens/'${ens}'/g'\
    | sed 's/fen/'${fen}'/g'\
    > getSubXdatafull_e$ens.py
done

# This section submits the python scripts on a HPC.
# Turned off in default
if [ 1 -eq 0 ];then
    rm -rf logs/*
    rm -rf submit_scripts/*e*.sh
    mkdir -p logs
    for ens in {${fen}..${len}}; do
        # Replace text in submit template file
        cat submit_scripts/submit_full.sh | sed 's/ens/'${ens}'/g' > submit_scripts/submit_e${ens}.sh
        bsub < submit_scripts/submit_e${ens}.sh
    done
fi
