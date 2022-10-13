# @Author: lshuns
# @Date:   2022-09-19 16:34:10
# @Last Modified by:   lshuns
# @Last Modified time: 2022-09-19 17:43:19

# usage: bias_estimate_script.py [-h] [--bias_type {alpha_LS,m_LS_pair,m_LS_tile}] [--in_file IN_FILE]
#                                [--out_path OUT_PATH] [--cols_e12 COLS_E12 COLS_E12]
#                                [--e_type {input,measured}] [--cols_g12 COLS_G12 COLS_G12]
#                                [--cols_PSFe12 COLS_PSFE12 COLS_PSFE12] [--col_weight COL_WEIGHT]
#                                [--col_label COL_LABEL] [--col_id_input COL_ID_INPUT]
#                                [--col_binning COL_BINNING]
#                                [--binning_edges [BINNING_EDGES [BINNING_EDGES ...]]]
#                                [--binning_min_max_Nbin BINNING_MIN_MAX_NBIN BINNING_MIN_MAX_NBIN BINNING_MIN_MAX_NBIN]
#                                [--binning_type {linear,log,quantile}] [--psf_frame]


# ###### v07D7: constShear whole (kids photometry)
# python bias_estimate_script.py --bias_type m_LS_tile\
#                                 --in_file /disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut.feather\
#                                 --out_path ./outputs/m_weiRaw_skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut.csv\
#                                 --cols_e12 e1_LF_r e2_LF_r\
#                                 --e_type measured\
#                                 --cols_g12 g1_in g2_in\
#                                 --col_weight oldweight_LF_r\
#                                 --col_label tile_label\
#                                 --col_binning Z_B\
#                                 --binning_edges 0.1 0.3 0.5 0.7 0.9 1.2 2.0
# binning on: Z_B
# provided bins [0.1 0.3 0.5 0.7 0.9 1.2 2. ]
# using tile based mCalFunc
# Number of sources in the catalogue 47870631
# selected objects (within bin edge) 47499698
# selected objects (weight>0) 38907653
# results saved to ./outputs/m_weiRaw_skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut.csv
# Elapsed:6:07.73,User=265.892,System=1526.923,CPU=487.5%.

###### v07D7: constShear whole (kids photometry), After AlphaRecal
python bias_estimate_script.py --bias_type m_LS_tile\
                                --in_file /disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut_A1_A2.feather\
                                --out_path ./outputs/m_skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut_A12.csv\
                                --cols_e12 AlphaRecalD2_e1 AlphaRecalD2_e2\
                                --e_type measured\
                                --cols_g12 g1_in g2_in\
                                --col_weight AlphaRecalC_weight\
                                --col_label tile_label\
                                --col_binning Z_B\
                                --binning_edges 0.1 0.3 0.5 0.7 0.9 1.2 2.0
# binning on: Z_B
# provided bins [0.1 0.3 0.5 0.7 0.9 1.2 2. ]
# using tile based mCalFunc
# Number of sources in the catalogue 38759673
# selected objects (within bin edge) 38759673
# selected objects (weight>0) 38759673
# results saved to ./outputs/m_skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut_A12.csv
# Elapsed:5:32.71,User=232.120,System=1327.859,CPU=468.8%.
