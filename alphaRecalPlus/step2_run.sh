# @Author: lshuns
# @Date:   2022-09-14 17:50:27
# @Last Modified by:   lshuns
# @Last Modified time: 2022-09-26 20:37:17

# ### running script
# usage: step2_methodD.py [-h] [--inpath INPATH] [--outDir OUTDIR]
#                         [--col_goldFlag COL_GOLDFLAG]
#                         [--col_weight COL_WEIGHT] [--col_snr COL_SNR]
#                         [--col_ZB COL_ZB] [--cols_e12 COLS_E12 COLS_E12]
#                         [--cols_psf_e12 COLS_PSF_E12 COLS_PSF_E12]
#                         [--Z_B_edges [Z_B_EDGES [Z_B_EDGES ...]]]

# step2_methodD.py: correct alpha in e1,2 with method D.

# optional arguments:
#   -h, --help            show this help message and exit
#   --inpath INPATH       the in path for the catalogue.
#                             supported formats: feather
#   --outDir OUTDIR       directory for the final catalogue.
#                             outpath name will be inpath_name + _A2 or _A2_goldclasses
#   --col_goldFlag COL_GOLDFLAG
#                         columns to the redshift gold class in the catalogue.
#   --col_weight COL_WEIGHT
#                         columns to the weight in the catalogue.
#   --col_snr COL_SNR     columns to the SNR in the catalogue.
#   --col_ZB COL_ZB       columns to the Z_B in the catalogue.
#   --cols_e12 COLS_E12 COLS_E12
#                         column names for e1_gal, e2_gal.
#   --cols_psf_e12 COLS_PSF_E12 COLS_PSF_E12
#                         column names for e1_psf, e2_psf.
#   --Z_B_edges [Z_B_EDGES [Z_B_EDGES ...]]
#                         edges for tomographic binning.

# >>>>>>>>>>>>>>>>> KiDS whole
python step2_methodD.py\
    --inpath /disks/shear16/ssli/KiDS/K1000_forSKiLLS/kids_photo_LF_321_shear_noSG_noWeiCut_newCut_817tiles_A1.feather\
    --outDir /disks/shear16/ssli/KiDS/K1000_forSKiLLS\
    --col_weight AlphaRecalC_weight\
    --col_snr model_SNratio\
    --col_ZB Z_B\
    --cols_e12 autocal_e1 autocal_e2\
    --cols_psf_e12 PSF_e1 PSF_e2\
    --Z_B_edges 0.1 0.3 0.5 0.7 0.9 1.2 2.0

# number original 31497505
# number after weight selection 26618226
# number after Z_B selection 26360327
# alpha map produced in 241.2777271270752 s
# number with meaningful e after D1 26360323 fraction 0.9999998482568141
# D1 finished in 143.0598111152649 s
# number with meaningful e after D2 26360323 fraction 1.0
# D2 finished in 704.8611278533936 s
# number in final cata 26360323
# final results saved to /disks/shear16/ssli/KiDS/K1000_forSKiLLS/kids_photo_LF_321_shear_noSG_noWeiCut_newCut_817tiles_A1_A2.feather
# Elapsed:20:40.25,User=1101.575,System=2720.411,CPU=308.1%.
