# @Author: lshuns
# @Date:   2023-02-06 11:26:31
# @Last Modified by:   lshuns
# @Last Modified time: 2023-03-08 10:15:35

# ### running script
# usage: step1_methodC.py [-h] [--inpath INPATH] [--outDir OUTDIR]
#                         [--col_weight COL_WEIGHT] [--col_var COL_VAR]
#                         [--col_snr COL_SNR] [--cols_e12 COLS_E12 COLS_E12]
#                         [--cols_e12_corr COLS_E12_CORR COLS_E12_CORR]
#                         [--cols_e12_raw COLS_E12_RAW COLS_E12_RAW]
#                         [--cols_psf_e12 COLS_PSF_E12 COLS_PSF_E12]

# step1_methodC.py: correct alpha in variance with method C.

# optional arguments:
#   -h, --help            show this help message and exit
#   --inpath INPATH       the in path for the catalogue.
#                             supported formats: feather
#   --outDir OUTDIR       directory for the final catalogue.
#                             outpath name will be inpath_name + _step1_methodC
#   --col_weight COL_WEIGHT
#                         columns to the weight in the catalogue.
#   --col_var COL_VAR     columns to the variance in the catalogue.
#   --col_snr COL_SNR     columns to the SNR in the catalogue.
#   --cols_e12 COLS_E12 COLS_E12
#                         column names for e1_gal, e2_gal.
#   --cols_e12_corr COLS_E12_CORR COLS_E12_CORR
#                         column names for e1, e2, correction.
#   --cols_e12_raw COLS_E12_RAW COLS_E12_RAW
#                         column names for e1, e2, raw measurement.
#                             this or cols_e12_corr should be exclusive.
#   --cols_psf_e12 COLS_PSF_E12 COLS_PSF_E12
#                         column names for e1_psf, e2_psf.

# # >>>>>>>>>>>>>>>>> SKiLLS 321 (KiDS photometry)
# python ../../alphaRecalPlus/step1_methodC.py\
#     --inpath /disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7p1_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut.feather\
#     --outDir /disks/shear16/ssli/ImSim/output/skills_v07D7\
#     --col_weight oldweight_LF_r\
#     --col_var LS_variance_LF_r\
#     --col_snr SNR_LF_r\
#     --cols_e12 e1_LF_r e2_LF_r\
#     --cols_e12_corr e1_corr_LF_r e2_corr_LF_r\
#     --cols_psf_e12 psf_e1_LF_r psf_e2_LF_r
# # number total 47869587
# # final results saved to /disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7p1_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut_A1.feather
# # Elapsed:30:56.56,User=1660.889,System=2302.149,CPU=213.4%.

# >>>>>>>>>>>>>>>>> SKiLLS 321 (Shark photometry)
python ../../alphaRecalPlus/step1_methodC.py\
    --inpath /disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7p1_LF_321_sharkPhotometry_shear_noSG_noWeiCut_newCut.feather\
    --outDir /disks/shear16/ssli/ImSim/output/skills_v07D7\
    --col_weight oldweight_LF_r\
    --col_var LS_variance_LF_r\
    --col_snr SNR_LF_r\
    --cols_e12 e1_LF_r e2_LF_r\
    --cols_e12_corr e1_corr_LF_r e2_corr_LF_r\
    --cols_psf_e12 psf_e1_LF_r psf_e2_LF_r
# number total 47928768
# sfinal results saved to /disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7p1_LF_321_sharkPhotometry_shear_noSG_noWeiCut_newCut_A1.feather
# Elapsed:32:29.74,User=1602.162,System=5014.516,CPU=339.3%.