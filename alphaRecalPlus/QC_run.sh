# @Author: lshuns
# @Date:   2022-10-20 16:02:25
# @Last Modified by:   lshuns
# @Last Modified time: 2022-10-20 22:40:22

#### for step1
# usage: QC_step1.py [-h] [--inpath INPATH] [--outpath OUTPATH] [--col_R COL_R]
#                    [--col_snr COL_SNR]
#                    [--cols_e12_raw COLS_E12_RAW COLS_E12_RAW]
#                    [--cols_psf_e12 COLS_PSF_E12 COLS_PSF_E12]
#                    [--col_var_before COL_VAR_BEFORE]
#                    [--col_var_after COL_VAR_AFTER]
#                    [--Nbins_R_snr NBINS_R_SNR NBINS_R_SNR]

# #### for step2
# usage: QC_step2.py [-h] [--inpath INPATH] [--outprefix OUTPREFIX]
#                    [--col_ZB COL_ZB] [--col_R COL_R] [--col_snr COL_SNR]
#                    [--cols_psf_e12 COLS_PSF_E12 COLS_PSF_E12]
#                    [--cols_e12_before COLS_E12_BEFORE COLS_E12_BEFORE]
#                    [--cols_e12_after COLS_E12_AFTER COLS_E12_AFTER]
#                    [--col_weight_after COL_WEIGHT_AFTER]
#                    [--Nbins_R_snr NBINS_R_SNR NBINS_R_SNR]
#                    [--Z_B_edges [Z_B_EDGES [Z_B_EDGES ...]]]

# >>>>>>>>>>>>>>>>> KiDS 321
#### step1
python QC_step1.py\
    --inpath /disks/shear16/ssli/KiDS/K1000_forSKiLLS/kids_photo_LF_321_shear_noSG_noWeiCut_newCut_817tiles_A1.feather\
    --outpath ./QC_alphaRecal_step1.png\
    --col_R R\
    --col_snr model_SNratio\
    --cols_e12_raw raw_e1 raw_e2\
    --cols_psf_e12 PSF_e1 PSF_e2\
    --col_var_before 2D_measurement_variance\
    --col_var_after AlphaRecalC_variance\
    --Nbins_R_snr 30 30
# QC plot for alphaRecal step 1 saved as ./QC_alphaRecal_step1.png
# Elapsed:3:36.00,User=650.283,System=3024.611,CPU=1701.3%.

#### step2
python QC_step2.py\
    --inpath /disks/shear16/ssli/KiDS/K1000_forSKiLLS/kids_photo_LF_321_shear_noSG_noWeiCut_newCut_817tiles_A1_A2.feather\
    --outprefix ./QC_alphaRecal_step2\
    --col_ZB Z_B\
    --col_R R\
    --col_snr model_SNratio\
    --cols_psf_e12 PSF_e1 PSF_e2\
    --cols_e12_before autocal_e1 autocal_e2\
    --cols_e12_after AlphaRecalD2_e1 AlphaRecalD2_e2\
    --col_weight_after AlphaRecalC_weight\
    --Nbins_R_snr 30 30\
    --Z_B_edges 0.1 0.3 0.5 0.7 0.9 1.2 2.0
# QC plot for alpha1 in SNR-R plane saved as ./QC_alphaRecal_step2_SNR_R_alpha1.png
# QC plot for alpha2 in SNR-R plane saved as ./QC_alphaRecal_step2_SNR_R_alpha2.png
# No handles with labels found to put in legend.
# QC plot for alpha in ZB bins saved as ./QC_alphaRecal_step2_ZBbins.png
# Elapsed:6:42.49,User=1277.586,System=6092.149,CPU=1831.0%.
