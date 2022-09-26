# @Author: lshuns
# @Date:   2022-05-31 15:07:11
# @Last Modified by:   lshuns
# @Last Modified time: 2022-09-26 18:46:02

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

# >>>>>>>>>>>>>>>>> KiDS 321
python step1_methodC.py\
    --inpath /disks/shear16/ssli/KiDS/K1000_forSKiLLS/kids_photo_LF_321_shear_noSG_noWeiCut_newCut_817tiles.feather\
    --outDir /disks/shear16/ssli/KiDS/K1000_forSKiLLS\
    --col_weight weight\
    --col_var 2D_measurement_variance\
    --col_snr model_SNratio\
    --cols_e12 autocal_e1 autocal_e2\
    --cols_e12_raw raw_e1 raw_e2\
    --cols_psf_e12 PSF_e1 PSF_e2
# number total 31497505
# final results saved to /disks/shear16/ssli/KiDS/K1000_forSKiLLS/kids_photo_LF_321_shear_noSG_noWeiCut_newCut_817tiles_A1.feather
# Elapsed:33:41.93,User=1206.974,System=2937.525,CPU=204.9%.

# >>>>>>>>>>>>>>>>> SKiLLS 321 (Shark photometry)
python step1_methodC.py\
    --inpath /disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7_LF_321_shear_noSG_noWeiCut_newCut.feather\
    --outDir /disks/shear16/ssli/ImSim/output/skills_v07D7\
    --col_weight oldweight_LF_r\
    --col_var LS_variance_LF_r\
    --col_snr SNR_LF_r\
    --cols_e12 e1_LF_r e2_LF_r\
    --cols_e12_corr e1_corr_LF_r e2_corr_LF_r\
    --cols_psf_e12 psf_e1_LF_r psf_e2_LF_r
# number total 47929683
# final results saved to /disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7_LF_321_shear_noSG_noWeiCut_newCut_A1.feather
# Elapsed:34:04.39,User=1694.625,System=4549.447,CPU=305.4%.

# >>>>>>>>>>>>>>>>> SKiLLS 321 (KiDS photometry)
python step1_methodC.py\
    --inpath /disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut.feather\
    --outDir /disks/shear16/ssli/ImSim/output/skills_v07D7\
    --col_weight oldweight_LF_r\
    --col_var LS_variance_LF_r\
    --col_snr SNR_LF_r\
    --cols_e12 e1_LF_r e2_LF_r\
    --cols_e12_corr e1_corr_LF_r e2_corr_LF_r\
    --cols_psf_e12 psf_e1_LF_r psf_e2_LF_r
# number total 47870631
# final results saved to /disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut_A1.feather
# Elapsed:31:52.78,User=1533.205,System=4291.438,CPU=304.5%.