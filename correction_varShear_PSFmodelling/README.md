# Correct shear biases from constant shear simulations

This repository contains codes for two corrections: one addressing the 'shear-interplay' effect and another for PSF modeling uncertainties.

## Account for the 'shear-interplay' effect 

[A_blending_fraction_2D.py](https://github.com/KiDS-WL/MultiBand_ImSim/tree/main/correction_varShear_PSFmodelling/A_blending_fraction_2D.py)

- Estimate the blending fraction in the 2D SNR and resolution bins for each tomographic bin. This requires the binning information defined during shear bias calculation from constant shear simulations (those from [bias_estimate_script_dataRewei.py](https://github.com/KiDS-WL/MultiBand_ImSim/tree/main/biasEstimation/bias_estimate_script_dataRewei.py)).

[B_dm_inBlendingSample_2D.py](https://github.com/KiDS-WL/MultiBand_ImSim/tree/main/correction_varShear_PSFmodelling/B_dm_inBlendingSample_2D.py)

- Calculate the shear bias difference between constant shear simulations and variable shear simulations in the 2D SNR and resolution bins for each tomographic bin using the blending-only samples.

[C_apply_dm_varShear_2D.py](https://github.com/KiDS-WL/MultiBand_ImSim/tree/main/correction_varShear_PSFmodelling/C_apply_dm_varShear_2D.py)

- Adjust the original shear biases to account for the dm estimated from variable shear simulations.

## Account for the PSF modelling uncertainties

[D_dm_PSFmodelling.py](https://github.com/KiDS-WL/MultiBand_ImSim/tree/main/correction_varShear_PSFmodelling/D_dm_PSFmodelling.py)

- Calculate the shear bias difference between idealised simulations and simulations with PSF modeling procedures.

[E_apply_dm_PSFmodelling.py](https://github.com/KiDS-WL/MultiBand_ImSim/tree/main/correction_varShear_PSFmodelling/E_apply_dm_PSFmodelling.py)

- Correct the shear biases to account for the dm estimated from the PSF modeling test.


