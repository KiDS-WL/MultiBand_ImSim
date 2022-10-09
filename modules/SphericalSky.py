# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-10-07 10:19:13
# @Last Modified by:   lshuns
# @Last Modified time: 2022-10-07 11:32:44

### functions related to the distance and area calculation in a sperical sky

import numpy as np

def RADECdistFunc(RADEC1, RADEC2, in_units = 'deg', out_units = 'deg'):
    """
    the distance between two points in (ra, dec)
    supported units: deg, arcmin, arcsec, radian
    """

    # avoid numerical issues
    RADEC1 = np.float64(RADEC1)
    RADEC2 = np.float64(RADEC2)

    # transfer everything to radian
    if in_units == 'deg':
        RA1 = RADEC1[0] * np.pi / 180. 
        DEC1 = RADEC1[1] * np.pi / 180. 
        RA2 = RADEC2[0] * np.pi / 180. 
        DEC2 = RADEC2[1] * np.pi / 180. 
    elif in_units == 'arcmin':
        RA1 = RADEC1[0] * np.pi / 60. / 180. 
        DEC1 = RADEC1[1] * np.pi / 60. / 180. 
        RA2 = RADEC2[0] * np.pi / 60. / 180. 
        DEC2 = RADEC2[1] * np.pi / 60. / 180. 
    elif in_units == 'arcsec':
        RA1 = RADEC1[0] * np.pi / 180. / 3600.
        DEC1 = RADEC1[1] * np.pi / 180. / 3600.
        RA2 = RADEC2[0] * np.pi / 180. / 3600.
        DEC2 = RADEC2[1] * np.pi / 180. / 3600.
    elif in_units == 'radian':
        RA1 = RADEC1[0]
        DEC1 = RADEC1[1] 
        RA2 = RADEC2[0] 
        DEC2 = RADEC2[1]
    else:
        raise Exception(f'Unsupported in_units: {in_units}')

    del RADEC1, RADEC2

    # calculate distance in radian
    ## using Eq (11) from https://www.aa.quae.nl/en/reken/afstanden.html
    tan2halfD =  (np.square(np.sin(0.5*(DEC2 - DEC1))) + np.cos(DEC1)*np.cos(DEC2)*np.square(np.sin(0.5*(RA2-RA1))))\
    /(np.square(np.sin(0.5*(DEC2 + DEC1))) + np.cos(DEC1)*np.cos(DEC2)*np.square(np.cos(0.5*(RA2-RA1))))
    del DEC2, DEC1, RA2, RA1
    ## back to distance
    Dfinal = 2 * np.arctan(np.sqrt(tan2halfD))
    del tan2halfD

    # transfer to outunits
    if out_units == 'deg':
        return Dfinal * 180. / np.pi 
    elif out_units == 'arcmin':
        return Dfinal * 180. * 60. / np.pi
    elif out_units == 'arcsec':
        return Dfinal * 180. * 3600. / np.pi
    elif out_units == 'radian':
        return Dfinal
    else:
        raise Exception(f'Unsupported out_units: {out_units}')

def RADECareaFunc(RA1, RA2, DEC1, DEC2, in_units='deg', out_units='deg'):
    """
    the sky area for a given range
    supported units: deg, arcmin, arcsec, radian
    """

    # avoid numerical issues
    RA1 = np.float64(RA1)
    RA2 = np.float64(RA2)
    DEC1 = np.float64(DEC1)
    DEC2 = np.float64(DEC2)

    # transfer everything to radian
    if in_units == 'deg':
        RA1 = RA1 * np.pi / 180. 
        DEC1 = DEC1 * np.pi / 180. 
        RA2 = RA2 * np.pi / 180. 
        DEC2 = DEC2 * np.pi / 180. 
    elif in_units == 'arcmin':
        RA1 = RA1 * np.pi / 60. / 180. 
        DEC1 = DEC1 * np.pi / 60. / 180. 
        RA2 = RA2 * np.pi / 60. / 180. 
        DEC2 = DEC2 * np.pi / 60. / 180. 
    elif in_units == 'arcsec':
        RA1 = RA1 * np.pi / 180. / 3600.
        DEC1 = DEC1 * np.pi / 180. / 3600.
        RA2 = RA2 * np.pi / 180. / 3600.
        DEC2 = DEC2 * np.pi / 180. / 3600.
    else:
        raise Exception(f'Unsupported in_units: {in_units}')

    # calculate area in radian
    area = (RA2 - RA1) * (np.sin(DEC2) - np.sin(DEC1))
    del DEC2, DEC1, RA2, RA1

    # transfer to outunits
    if out_units == 'deg':
        return area * 180.**2 / np.pi ** 2
    elif out_units == 'arcmin':
        return area * 180.**2 * 60.**2 / np.pi**2
    elif out_units == 'arcsec':
        return area * 180.**2 * 3600.**2 / np.pi**2
    elif out_units == 'radian':
        return area
    else:
        raise Exception(f'Unsupported out_units: {out_units}')

