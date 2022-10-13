# alphaRecalPlus

Empirically calibrate the PSF contaminations in the *lens*fit measurements (weight and ellipticity).

[step1_methodC.py](https://github.com/KiDS-WL/MultiBand_ImSim/blob/main/alphaRecalPlus/step1_methodC.py)

- The first part code calibrating the weight

- It should be applied to the combined (mosaic) catalogue after photometric and *lens*fit selections except for the weight cut.

[step2_methodD.py](https://github.com/KiDS-WL/MultiBand_ImSim/blob/main/alphaRecalPlus/step2_methodD.py)

- The second part code calibrating the measured ellipticity

- It should be applied after any redshift calibrations and tomographic binning.