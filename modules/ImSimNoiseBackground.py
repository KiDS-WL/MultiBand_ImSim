# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2020-11-26 15:44:09
# @Last Modified by:   lshuns
# @Last Modified time: 2021-01-29 16:23:15

### Everything about the background noise

import galsim
import numpy as np

def GaussianNoise(noise_rms, rng_seed=940120):
    """
    Generate Gaussian background noise.

    Parameters
    ----------
    noise_rms : float
        Observational noise rms for the whole image.   
    rng_seed : int, optional (default: 940120) 
        Seed for random number generator.

    Returns
    -------
    noise: galsim object
        Background noise.
    """

    rng = galsim.BaseDeviate(rng_seed)
    noise = galsim.GaussianNoise(rng, sigma=noise_rms)

    return noise