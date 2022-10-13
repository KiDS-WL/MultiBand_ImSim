# Post-processing of the MultiBand_ImSim pipeline outputs

[A_combine_everything.py](https://github.com/KiDS-WL/MultiBand_ImSim/blob/main/post_processing/skills_v07D7/A_combine_everything.py)

- To combine everything into one big catalogue.

[B_skills_flag.py](https://github.com/KiDS-WL/MultiBand_ImSim/blob/main/post_processing/skills_v07D7/B_skills_flag.py)

- To add flag columns for KiDS-like photo and lensfit selections.

- NOTE: it only adds flag columns but does not remove anything.

[C_split_shear_rot.py](https://github.com/KiDS-WL/MultiBand_ImSim/blob/main/post_processing/skills_v07D7/C_split_shear_rot.py)

- To split the combined catalogue into subsamples, each with unique inputs.