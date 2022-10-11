# MultiBand_ImSim

The Multi-Band Image Simulation pipeline for generating multi-band images and creating joint redshift-shear mock catalogues.

It contains 8 modules labelled as

- taskID 1: simulate images
- taskID 2: swarp images
- taskID 3: detect objects
- taskID 4: measure photometry
- taskID 5: measure photo-z
- taskID 6_1: psf modelling
- taskID 6_2: measure shapes
- taskID 7: create a combined catalogue

For detailed usage information, run

    `python modules/Run.py --help`

For some test scripts, look into [scripts](https://github.com/KiDS-WL/MultiBand_ImSim/tree/main/scripts)

For the correction of PSF contamination in *lens*fit measurements, look into [alphaRecalPlus](https://github.com/KiDS-WL/MultiBand_ImSim/tree/main/alphaRecalPlus)

## Link to the input catalogues

* [SURFS-based galaxies](https://surfdrive.surf.nl/files/index.php/s/uegK5tc15TbWib7)
* [Trilegal stars](https://surfdrive.surf.nl/files/index.php/s/dMsqnkeEUFSSLHE)

## Link to the SKiLLS output catalogues

SKiLLS: SURFS-based KiDS-Legacy-Like Simulations

### SKiLLS nine-band outputs

For the generation code, look into [post_processing/skills_v07D7](https://github.com/KiDS-WL/MultiBand_ImSim/tree/main/post_processing/skills_v07D7)

For catalogues, download [here](https://surfdrive.surf.nl/files/index.php/s/iSvDmHQJjDa0ewG)

> NOTE: these catalogues contain all objects detected by SExtractor, including false detections and stars.
>
> The catalogues contain flag columns for applying KiDS-like photometry / *lens*fit selections (see readme.txt within the shared folder for details)
