#!/bin/bash

#remove all copies of extra files
#ETO and soil shortwave radiation (dswrf)
cd /home/kdl/Insync/OneDrive/NRT_CPC_Internship/Data/SubX/GMAO

for f in ETo*.nc4;do
    if [[  "${#f}" -gt 19  ]];then
    rm "${f}"; fi
done

for f in dswrf*.nc;do
    if [[  "${#f}" -gt 24  ]];then
    rm "${f}"; fi
done

#EDDI
cd /home/kdl/Insync/OneDrive/NRT_CPC_Internship/Data/EDDI/convert_2_nc

for f in EDDI_ET*.nc4;do
    if [[  "${#f}" -gt 28  ]];then
    rm "${f}"; fi
done

for f in eddi_mer*.nc4;do
    if [[  "${#f}" -gt 16  ]];then
    rm "${f}"; fi
done

#SMERGE soil moisture
cd /home/kdl/Insync/OneDrive/NRT_CPC_Internship/Data/SMERGE_SM/Raw_data

for f in 1*.nc4;do
    if [[  "${#f}" -gt 13  ]];then
    rm "${f}"; fi
done

for f in 2*.nc4;do
    if [[  "${#f}" -gt 13  ]];then
    rm "${f}"; fi
done

#EDDI.asc files
cd /home/kdl/Insync/OneDrive/NRT_CPC_Internship/Data/EDDI

for f in *.asc;do
    if [[  "${#f}" -gt 28  ]];then
    rm "${f}"; fi
done


#SMERGE anomalies converted to SubX format
cd /home/kdl/Insync/OneDrive/NRT_CPC_Internship/Data/SMERGE_SM/SM_SubX_values

for f in *.nc4;do
    if [[  "${#f}" -gt 31  ]];then
    rm "${f}"; fi
done

#return to script directory
cd /home/kdl/Insync/OneDrive/NRT_CPC_Internship/Scripts
