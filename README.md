# This code was for a research project entitled "Subseasonal forecast skill of evaporative demand, soil moisture, and flash drought onset from two dynamic models over the contiguous United States". 

## There are scripts which are for downloading data and located in directory [download_data] (https://github.com/kdl0013/Flash-Drought-Forecast-Skill-GEFSv12-and-GMAOv2p1/tree/main/download_data)
0AA_download_SubX_models.sh
0AB_download_observations.sh
0A_download_data.sh
0A_download_data_GEFSv12.sh
0b_download_metdata_wget.sh

wget_MERRA.sh
wget_GMAO_2021.sh
wget_MERRA_small.sh

00_get_Priestley_Taylor_coefficients.py
00_download_ERA5_Land.py
00_download_MERRA2.py
00_download_SubX.py
00_create_wget_GEOS_2017_2021_years.py
00_merge_GMAO_models_2017_2021.py

01_combine_models_GEFSv12.py
01a_day_average_GEFSv12_HPC.py
01a_day_average_GEFSv12_soil_moisture_only.py
01c_remove_S_dimension_GEFSv12.py
01c_restrict_to_CONUS_GEFSv12.py
01d_remove_S_dimension_GEFSv12.py
01d_select_RSMAS_dates.py
01e_remove_S_dimsension_all_models.py


0a_download_SMERGE_RZSM.py
0c_download_EDDI.py
0c_download_EDDI_wget.sh
0d_download_EDDI_wget.sh


## Then there are scripts which further process data and located in [process_data] (https://github.com/kdl0013/Flash-Drought-Forecast-Skill-GEFSv12-and-GMAOv2p1/tree/main/process_data)
1A_preprocess_data.sh
1b1_make_reference_ET.py
1b2_convert_RZSM_to_m3_m3.py
1b3_remove_S_dimsension.py
1b4_make_empty_anomaly_files.py
1b_ETo_and_make_empty_netcdf_files.py
1c2_anomaly_RZSM_mean_for_ACC.py
1c2_make_anomaly_MME.py

1c_anomaly_RZSM.py
1c_create_mean_EMC.py
1c_make_anomaly.py
1c_make_anomaly_EMC.py
1c_make_anomaly_GEFSv12.py
1c_make_anomaly_GMAO.py

1d_EDDI.py
1d_anomaly_ETo.py
1e_EDDI.py
1e_anomaly_RZSM.py
1e_check_for_missing_data.py
1f_anomaly_ETo.py
1f_fix_individual_anomalies.py





