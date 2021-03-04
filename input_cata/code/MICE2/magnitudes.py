# @Author: ssli
# @Date:   2021-02-25, 13:10:58
# @Last modified by:   ssli
# @Last modified time: 2021-02-25, 17:26:53

### apply magnitude correction to mice2

import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.table import Table

def mag_correction(mag, redshift):
    """
    Evolution correction as described in the official MICE2 manual
    https://www.dropbox.com/s/0ffa8e7463n8h1q/README_MICECAT_v2.0_for_new_CosmoHub.pdf?dl=0

    Parameters
    ----------
    mag : array_like
        Uncorrected model magnitudes.
    redshift : array_like
        True galaxy redshift (z_cgal).

    Returns
    -------
    mag_evo : array_like
        Evolution corrected model magnitudes.
    """
    return mag - 0.8 * (np.arctan(1.5 * redshift) - 0.1489)


def magnification_correction(mag, kappa):
    """
    Magnification calculated from the convergence, following Fosalba+15 eq. 21.

    Parameters
    ----------
    mag : array_like
        (Evolution corrected) model magnitudes.
    kappa : array_like
        Convergence field at the galaxy positions.

    Returns
    -------
    mag_magnified : array_like
        Magnitudes corrected for magnification.
    """
    d_mu = 2.0 * kappa
    mag_magnified = mag - 2.5 * np.log10(1 + d_mu)

    return mag_magnified


if __name__=='__main__':

    #################### test catalogue (4deg)
    ## output
    out_file = '/disks/shear10/ssli/ImSim/input/MICE2_cata/MICE2_ra35-37_dec5-7_corrected.fits'

    ## input
    in_file = '/disks/shear10/ssli/ImSim/input/MICE2_cata/MICE2_ra35-37_dec5-7.fits'
    with fits.open(in_file) as hdul:
        MICE_cata = Table(hdul[1].data)

    ## useful values
    kappa = MICE_cata['kappa']
    redshift = MICE_cata['z_cgal']

    ## bands
    filters = ['sdss_u_true', 'sdss_g_true', 'sdss_r_true', 'sdss_i_true', 'sdss_z_true',
                'des_asahi_full_y_true', 'vhs_j_true', 'vhs_h_true', 'vhs_ks_true',
                'cosmos_subaru_r_true', 'cosmos_subaru_i_true']
    abs_filters = ['sdss_r_abs_mag', 'cosmos_subaru_r_abs_mag', 'cosmos_subaru_i_abs_mag']

    ## for apparent magnitude
    for filter in filters:
        mag_ori = MICE_cata[filter]

        # evolution correction
        mag_evo = mag_correction(mag_ori, redshift)
        # magnification correction
        mag_lensed = magnification_correction(mag_evo, kappa)
        MICE_cata.add_column(mag_lensed, name=f'{filter}_lensed')
        MICE_cata.remove_column(filter)

    ## for absolute magnitude
    for filter in abs_filters:
        mag_ori = MICE_cata[filter]

        # evolution correction
        mag_evo = mag_correction(mag_ori, redshift)
        MICE_cata.add_column(mag_evo, name=f'{filter}_evo')

        # magnification correction
        mag_lensed = magnification_correction(mag_evo, kappa)
        MICE_cata.add_column(mag_lensed, name=f'{filter}_lensed')


    MICE_cata.write(out_file)
