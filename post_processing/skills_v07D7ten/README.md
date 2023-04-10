# Post-processing for SKiLLS v07D7ten catalogues

This repository contains codes for conducting post-processing on the ten-band SKiLLS catalogues.

[A_combine_everything.py](https://github.com/KiDS-WL/MultiBand_ImSim/blob/main/post_processing/skills_v07D7ten/A_combine_everything.py)

- Combine all data into a single, comprehensive catalogue.

[B_skills_flag.py](https://github.com/KiDS-WL/MultiBand_ImSim/blob/main/post_processing/skills_v07D7ten/B_skills_flag.py)

- Add flag columns for KiDS-like photometry and *lens*fit selections.

- NOTE: This step only adds flag columns without removing any data.

[C_apply_shear_flag.py](https://github.com/KiDS-WL/MultiBand_ImSim/blob/main/post_processing/skills_v07D7ten/C_apply_shear_flag.py)

- Apply *lens*fit-shear-catalogue selections (excluding the weight cut) in preparation for alphaRecal.

[D_alphaRecal_step1.sh](https://github.com/KiDS-WL/MultiBand_ImSim/blob/main/post_processing/skills_v07D7ten/D_alphaRecal_step1.sh)

- Execute the first step of AlphaRecal to obtain the corrected LF weights (AlphaRecalC_weight).

[E_combineAlphaRecal_split_shear_rot.py](https://github.com/KiDS-WL/MultiBand_ImSim/blob/main/post_processing/skills_v07D7ten/E_combineAlphaRecal_split_shear_rot.py)

- Assign AlphaRecal results back to the entire catalogue and divide it into subsamples, each with unique input parameters.

> The catalogues can be accessed through this [link](https://surfdrive.surf.nl/files/index.php/s/iSvDmHQJjDa0ewG?path=%2Fskills_v07D7ten)