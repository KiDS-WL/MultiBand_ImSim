# @Author: lshuns
# @Date:   2023-02-04 16:04:20
# @Last Modified by:   lshuns
# @Last Modified time: 2023-04-11 14:08:00

# usage: bias_estimate_script.py [-h] [--bias_type {alpha_LS,m_LS_pair,m_LS_tile}] [--in_file IN_FILE]
#                                [--out_path OUT_PATH] [--cols_e12 COLS_E12 COLS_E12]
#                                [--e_type {input,measured}] [--cols_g12 COLS_G12 COLS_G12]
#                                [--cols_PSFe12 COLS_PSFE12 COLS_PSFE12] [--col_weight COL_WEIGHT]
#                                [--col_label COL_LABEL] [--col_id_input COL_ID_INPUT]
#                                [--col_binning COL_BINNING]
#                                [--binning_edges [BINNING_EDGES [BINNING_EDGES ...]]]
#                                [--binning_min_max_Nbin BINNING_MIN_MAX_NBIN BINNING_MIN_MAX_NBIN BINNING_MIN_MAX_NBIN]
#                                [--binning_type {linear,log,quantile}] [--psf_frame]

# ###### v07D7p1: constShear whole (kids photometry)
python bias_estimate_script.py --bias_type m_LS_tile\
                                --in_file /disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7p1_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut.feather\
                                --out_path ./results/m_weiRaw_skills_v07D7p1_kidsPhotometry.csv\
                                --cols_e12 e1_LF_r e2_LF_r\
                                --e_type measured\
                                --cols_g12 g1_in g2_in\
                                --col_weight oldweight_LF_r\
                                --col_label tile_label\
                                --col_binning Z_B\
                                --binning_edges 0.1 0.3 0.5 0.7 0.9 1.2 2.0
# # ### out info
# binning on: Z_B
# provided bins [0.1 0.3 0.5 0.7 0.9 1.2 2. ]
# using tile based mCalFunc
# Number of sources in the catalogue 47869587
# selected objects (within bin edge) 47510465
# selected objects (weight>0) 38916184
# results saved to ./results/m_weiRaw_skills_v07D7p1_kidsPhotometry.csv
# Elapsed:5:41.48,User=270.147,System=1405.989,CPU=490.8%.