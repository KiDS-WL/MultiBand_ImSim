# Estimate shear bias using surface of doom

This repository provides scripts to estimate shear biases for any given sample using the surface of doom. 

The process involves two steps:

- Generating the surface of doom, which is done using the scripts with prefix from A to B. However, if the required surface of doom is already present in the [SurfaceOfDoom](https://github.com/KiDS-WL/MultiBand_ImSim/tree/main/biasEstimationWithSurfaceOfDoom/SurfaceOfDoom) folder, there is no need to run these scripts.

- Calculating *m* for any given sample using the function provided in [m_from_doom_func.py](https://github.com/KiDS-WL/MultiBand_ImSim/tree/main/biasEstimationWithSurfaceOfDoom/m_from_doom_func.py). Scripts with prefix C show some examples of how to use this function.