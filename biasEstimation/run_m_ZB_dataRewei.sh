# @Author: lshuns
# @Date:   2022-08-17 15:04:38
# @Last Modified by:   lshuns
# @Last Modified time: 2022-09-29 13:00:25

# usage: bias_estimate_script_dataRewei.py [-h]
#                                          [--bias_type {m_LS_DataRewei_2D}]
#                                          [--in_file IN_FILE]
#                                          [--out_path OUT_PATH]
#                                          [--cols_e12 COLS_E12 COLS_E12]
#                                          [--cols_g12 COLS_G12 COLS_G12]
#                                          [--col_weight_sim_data COL_WEIGHT_SIM_DATA COL_WEIGHT_SIM_DATA]
#                                          [--col_label COL_LABEL]
#                                          [--col_id_input COL_ID_INPUT]
#                                          [--col_binning_sim_data COL_BINNING_SIM_DATA COL_BINNING_SIM_DATA]
#                                          [--binning_edges [BINNING_EDGES [BINNING_EDGES ...]]]
#                                          [--binning_min_max_Nbin BINNING_MIN_MAX_NBIN BINNING_MIN_MAX_NBIN BINNING_MIN_MAX_NBIN]
#                                          [--binning_type {linear,log,quantile}]
#                                          [--psf_frame]
#                                          [--in_file_data IN_FILE_DATA]
#                                          [--bin1_info BIN1_INFO BIN1_INFO BIN1_INFO]
#                                          [--bin2_info BIN2_INFO BIN2_INFO BIN2_INFO]
#                                          [--fitting_method {tile_based,pair_based}]
#                                          [--save_surface_prefix SAVE_SURFACE_PREFIX]
#                                          [--save_bounds_prefix SAVE_BOUNDS_PREFIX]

###### v07D7: constShear whole (kids photometry)
python bias_estimate_script_dataRewei.py --bias_type m_LS_DataRewei_2D\
                                --in_file /disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut.feather\
                                --out_path ./outputs/m_weiRaw_DR4817Rewei_skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut.csv\
                                --cols_e12 e1_LF_r e2_LF_r\
                                --cols_g12 g1_in g2_in\
                                --col_weight_sim_data oldweight_LF_r weight\
                                --col_label tile_label\
                                --col_binning_sim_data Z_B Z_B\
                                --binning_edges 0.1 0.3 0.5 0.7 0.9 1.2 2.0\
                                --in_file_data /disks/shear16/ssli/KiDS/K1000_forSKiLLS/kids_photo_LF_321_noSG_noWeiCut_newCut_817tiles.feather\
                                --bin1_info SNR_LF_r model_SNratio 20\
                                --bin2_info R R 20\
                                --fitting_method tile_based\
                                --save_surface_prefix ./outputs/m_weiRaw_DR4817Rewei_skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut\
                                --save_bounds_prefix ./outputs/m_weiRaw_DR4817Rewei_skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut
# ### out info
# binning on: ['Z_B', 'Z_B']
# provided bins [0.1 0.3 0.5 0.7 0.9 1.2 2. ]
# using data reweighting in 2 parameters
#      with the tile based mCalFunc
# Number of sources in the data 31497513
#      within the edges 31200045
#      with weight > 0 26453911
# Number of sources in the simulation 47870631
#      within the edges 47499698
#      with weight > 0 38907653
# bin1_bounds saved to ./outputs/m_weiRaw_DR4817Rewei_skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut_whole_bin1.npy
# bin1_bounds saved to ./outputs/m_weiRaw_DR4817Rewei_skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut_whole_bin2.npy
# the bias surface saved as ./outputs/m_weiRaw_DR4817Rewei_skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut_whole.csv
# bin1_bounds saved to ./outputs/m_weiRaw_DR4817Rewei_skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut_bin0_bin1.npy
# bin1_bounds saved to ./outputs/m_weiRaw_DR4817Rewei_skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut_bin0_bin2.npy
# the bias surface saved as ./outputs/m_weiRaw_DR4817Rewei_skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut_bin0.csv
# bin1_bounds saved to ./outputs/m_weiRaw_DR4817Rewei_skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut_bin1_bin1.npy
# bin1_bounds saved to ./outputs/m_weiRaw_DR4817Rewei_skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut_bin1_bin2.npy
# the bias surface saved as ./outputs/m_weiRaw_DR4817Rewei_skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut_bin1.csv
# bin1_bounds saved to ./outputs/m_weiRaw_DR4817Rewei_skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut_bin2_bin1.npy
# bin1_bounds saved to ./outputs/m_weiRaw_DR4817Rewei_skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut_bin2_bin2.npy
# the bias surface saved as ./outputs/m_weiRaw_DR4817Rewei_skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut_bin2.csv
# bin1_bounds saved to ./outputs/m_weiRaw_DR4817Rewei_skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut_bin3_bin1.npy
# bin1_bounds saved to ./outputs/m_weiRaw_DR4817Rewei_skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut_bin3_bin2.npy
# the bias surface saved as ./outputs/m_weiRaw_DR4817Rewei_skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut_bin3.csv
# bin1_bounds saved to ./outputs/m_weiRaw_DR4817Rewei_skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut_bin4_bin1.npy
# bin1_bounds saved to ./outputs/m_weiRaw_DR4817Rewei_skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut_bin4_bin2.npy
# the bias surface saved as ./outputs/m_weiRaw_DR4817Rewei_skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut_bin4.csv
# bin1_bounds saved to ./outputs/m_weiRaw_DR4817Rewei_skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut_bin5_bin1.npy
# bin1_bounds saved to ./outputs/m_weiRaw_DR4817Rewei_skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut_bin5_bin2.npy
# the bias surface saved as ./outputs/m_weiRaw_DR4817Rewei_skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut_bin5.csv
# results saved to ./outputs/m_weiRaw_DR4817Rewei_skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut.csv
# Elapsed:12:38.16,User=593.538,System=1898.794,CPU=328.7%.

