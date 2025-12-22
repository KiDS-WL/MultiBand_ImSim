# MultiBand_ImSim

A pipeline designed to generate multi-band images and create joint redshift-shear mock catalogues.

It contains 8 modules labelled as

- taskID 1: simulate images
- taskID 2: swarp images
- taskID 3: detect objects
- taskID 4: measure photometry
- taskID 5: estimate photo-z
- taskID 6_1: model PSF
- taskID 6_2: measure shapes
- taskID 7: create a combined catalogue

For usage instructions, run

    `python modules/Run.py --help`

For sample test scripts, refer to the [scripts](https://github.com/KiDS-WL/MultiBand_ImSim/tree/main/scripts) directory.

## Link to the input catalogues

* [SURFS-based galaxies](https://drive.google.com/file/d/1lxaMGsKBWEScJfluMvhZAquceK6ln5u0/view?usp=sharing)
* [Trilegal stars](https://drive.google.com/file/d/1cH-0orCV-uASEkThbRwTvMJHvOK5Xt3a/view?usp=sharing)

## Link to the SKiLLS output catalogues

SKiLLS: SURFS-based KiDS-Legacy-Like Simulations

The post-processing codes can be found in the [post_processing](https://github.com/KiDS-WL/MultiBand_ImSim/tree/main/post_processing) directory.

The catalogues can be download [here](https://drive.google.com/drive/folders/1IkhFhU73mtNoJ6bMU8XtidrxVhrzsmNn?usp=sharing)

> NOTE: These catalogues include all objects detected by SExtractor, encompassing false detections and stars. However, they also feature flag columns for applying KiDS-like photometry and *lens*fit selections. For more details, please refer to the readme.txt file within the shared folder.

## AlphaRecal for *lens*fit measurements

The raw measurements from *lens*fit exhibit biases, primarily resulting from PSF anisotropy. We address these PSF-related distortions using an empirical approach, as elaborated in Section 4 of [Li et al. (2023)](https://ui.adsabs.harvard.edu/abs/2023A%26A...670A.100L/abstract).

The relevant codes can be found in the [alphaRecalPlus](https://github.com/KiDS-WL/MultiBand_ImSim/tree/main/alphaRecalPlus) repository.

## SDSS filters to KiDS/VIKING filters

After generating the images, we discovered that our SURFS-based galaxies possess photometry derived from SDSS filters, which differ slightly from the KiDS/VIKING filters. Consequently, we applied an empirical correction to the photometry measured from the simulated images to account for the filter discrepancies.

A comprehensive explanation is available in Appendix C of [Li et al. (2023)](https://ui.adsabs.harvard.edu/abs/2023A%26A...670A.100L/abstract), and the pertinent codes can be found in the [change_filters](https://github.com/KiDS-WL/MultiBand_ImSim/tree/main/change_filters) repository.

## Shear calibration 

Our shear calibration contains two main steps: 

1. Estimate shear biases using constant shear image simulations, adhering to previous KiDS conventions, and apply data reweighting based on the lensfit-reported model signal-to-noise ratio and object resolution to accommodate potential discrepancies between simulations and data. Further information is provided in Section 5.1 of [Li et al. (2023)](https://ui.adsabs.harvard.edu/abs/2023A%26A...670A.100L/abstract), with the associated codes located in the [biasEstimation](https://github.com/KiDS-WL/MultiBand_ImSim/tree/main/biasEstimation) repository.

2. Refine the shear biases from constant shear image simulations to address high-order effects stemming from the 'shear-interplay' effect and PSF modeling uncertainties. Further information can be found in Sections 5.2 and 5.3 of [Li et al. (2023)](https://ui.adsabs.harvard.edu/abs/2023A%26A...670A.100L/abstract), with the corresponding codes located in the [correction_varShear_PSFmodelling](https://github.com/KiDS-WL/MultiBand_ImSim/tree/main/correction_varShear_PSFmodelling) repository.

## Reference

Should the code or catalogues prove beneficial to your projects, please consider citing our paper: [Li et al. 2023](https://ui.adsabs.harvard.edu/abs/2023A%26A...670A.100L/abstract).