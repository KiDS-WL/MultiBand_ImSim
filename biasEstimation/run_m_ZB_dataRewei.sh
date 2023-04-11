# @Author: lshuns
# @Date:   2023-02-04 16:04:35
# @Last Modified by:   lshuns
# @Last Modified time: 2023-04-11 13:53:47

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

###### v07D7p1: constShear whole (kids photometry)
python bias_estimate_script_dataRewei.py --bias_type m_LS_DataRewei_2D\
                                --in_file /disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7p1_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut.feather\
                                --out_path ./results/m_weiRaw_K1000_LF_321_skills_v07D7p1_kidsPhotometry.csv\
                                --cols_e12 e1_LF_r e2_LF_r\
                                --cols_g12 g1_in g2_in\
                                --col_weight_sim_data oldweight_LF_r weight\
                                --col_label tile_label\
                                --col_binning_sim_data Z_B Z_B\
                                --binning_edges 0.1 0.3 0.5 0.7 0.9 1.2 2.0\
                                --in_file_data /disks/shear10/ssli/KiDS/K1000_LF_321_mosaic/KIDS_1006tiles_r_SDSS.V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_selec_noWeiCut.feather\
                                --bin1_info SNR_LF_r model_SNratio 20\
                                --bin2_info R R 20\
                                --fitting_method tile_based\
                                --save_surface_prefix ./results/m_weiRaw_K1000_LF_321_skills_v07D7p1_kidsPhotometry\
                                --save_bounds_prefix ./results/m_weiRaw_K1000_LF_321_skills_v07D7p1_kidsPhotometry

# # ### out info
# binning on: ['Z_B', 'Z_B']
# provided bins [0.1 0.3 0.5 0.7 0.9 1.2 2. ]
# using data reweighting in 2 parameters
#      with the tile based mCalFunc
# Number of sources in the data 38026162
#      within the edges 37667985
#      with weight > 0 31869107
# Number of sources in the simulation 47869587
#      within the edges 47510465
#      with weight > 0 38916184
# bin1_bounds saved to ./results/m_weiRaw_K1000_LF_321_skills_v07D7p1_kidsPhotometry_whole_bin1.npy
# bin1_bounds saved to ./results/m_weiRaw_K1000_LF_321_skills_v07D7p1_kidsPhotometry_whole_bin2.npy
# the bias surface saved as ./results/m_weiRaw_K1000_LF_321_skills_v07D7p1_kidsPhotometry_whole.csv
# bin1_bounds saved to ./results/m_weiRaw_K1000_LF_321_skills_v07D7p1_kidsPhotometry_bin0_bin1.npy
# bin1_bounds saved to ./results/m_weiRaw_K1000_LF_321_skills_v07D7p1_kidsPhotometry_bin0_bin2.npy
# the bias surface saved as ./results/m_weiRaw_K1000_LF_321_skills_v07D7p1_kidsPhotometry_bin0.csv
# bin1_bounds saved to ./results/m_weiRaw_K1000_LF_321_skills_v07D7p1_kidsPhotometry_bin1_bin1.npy
# bin1_bounds saved to ./results/m_weiRaw_K1000_LF_321_skills_v07D7p1_kidsPhotometry_bin1_bin2.npy
# the bias surface saved as ./results/m_weiRaw_K1000_LF_321_skills_v07D7p1_kidsPhotometry_bin1.csv
# bin1_bounds saved to ./results/m_weiRaw_K1000_LF_321_skills_v07D7p1_kidsPhotometry_bin2_bin1.npy
# bin1_bounds saved to ./results/m_weiRaw_K1000_LF_321_skills_v07D7p1_kidsPhotometry_bin2_bin2.npy
# the bias surface saved as ./results/m_weiRaw_K1000_LF_321_skills_v07D7p1_kidsPhotometry_bin2.csv
# bin1_bounds saved to ./results/m_weiRaw_K1000_LF_321_skills_v07D7p1_kidsPhotometry_bin3_bin1.npy
# bin1_bounds saved to ./results/m_weiRaw_K1000_LF_321_skills_v07D7p1_kidsPhotometry_bin3_bin2.npy
# the bias surface saved as ./results/m_weiRaw_K1000_LF_321_skills_v07D7p1_kidsPhotometry_bin3.csv
# bin1_bounds saved to ./results/m_weiRaw_K1000_LF_321_skills_v07D7p1_kidsPhotometry_bin4_bin1.npy
# bin1_bounds saved to ./results/m_weiRaw_K1000_LF_321_skills_v07D7p1_kidsPhotometry_bin4_bin2.npy
# the bias surface saved as ./results/m_weiRaw_K1000_LF_321_skills_v07D7p1_kidsPhotometry_bin4.csv
# bin1_bounds saved to ./results/m_weiRaw_K1000_LF_321_skills_v07D7p1_kidsPhotometry_bin5_bin1.npy
# bin1_bounds saved to ./results/m_weiRaw_K1000_LF_321_skills_v07D7p1_kidsPhotometry_bin5_bin2.npy
# the bias surface saved as ./results/m_weiRaw_K1000_LF_321_skills_v07D7p1_kidsPhotometry_bin5.csv
# results saved to ./results/m_weiRaw_K1000_LF_321_skills_v07D7p1_kidsPhotometry.csv
# Elapsed:16:33.36,User=646.704,System=2320.450,CPU=298.6%.