# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-08-30 17:02:54
# @Last Modified by:   lshuns
# @Last Modified time: 2022-09-29 20:07:40

### compare ZB results between SKiLLS (transform before and after) and KiDS

import pandas as pd 
import numpy as np

import plotting

# ++++++++++++++++++++++++++ I/O

outpath = 'show'
outpath = './plots/ZBcomp.pdf'

# data
infile = '/disks/shear16/ssli/KiDS/K1000_forSKiLLS/kids_photo_LF_321_photo_noSG_817tiles.feather'
data_whole = pd.read_feather(infile)
## get galaxy only
data_whole = data_whole[(data_whole['SG_FLAG']==1) & (data_whole['SG2DPHOT']==0)]
data_whole.reset_index(drop=True, inplace=True)
## used columns
data_whole = data_whole[['Z_B']]
## sky area
area_data = 636.2882737918676

# simulation before transform
infile = '/disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7_LF_321_photo_noSG_shear_m283m283_rot_0.feather'
sim_whole0 = pd.read_feather(infile)
## get galaxy only
sim_whole0 = sim_whole0[sim_whole0['id_input']>0]
sim_whole0.reset_index(drop=True, inplace=True)
## used columns
sim_whole0 = sim_whole0[['Z_B', 'MAG_AUTO']]
## sky area
area_sim = 108.

# simulation after transform
infile = '/disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7_LF_321_kidsPhotometry_photo_noSG_shear_m283m283_rot_0.feather'
sim_whole = pd.read_feather(infile)
## get galaxy only
sim_whole = sim_whole[sim_whole['id_input']>0]
sim_whole.reset_index(drop=True, inplace=True)
## used columns
sim_whole = sim_whole[['Z_B', 'MAG_AUTO']]
## sky area
area_sim = 108.

########## general 
font_size = 16

# ++++++++++++++++++++++++++ plot
paras = [data_whole['Z_B'].values, sim_whole['Z_B'].values, sim_whole0['Z_B'].values]
wgs = [np.ones(len(data_whole))/area_data, np.ones(len(sim_whole))/area_sim, np.ones(len(sim_whole0))/area_sim]
COLORs = ['r', 'b', 'green']
LABELs = ['KiDS DR4', 'SKiLLS (KiDS)', 'SKiLLS (SDSS)']
LINEs = ['-', '-', '--']
LINEWs = [2, 2, 2]

XRANGE = [0., 2.]
dz = 0.1
nbins = int((XRANGE[1]-XRANGE[0])/dz)
xtick_spe = [[0.5, 1.0, 1.5, 2.0], ['0.5', '1.0', '1.5', '2.0']]

ytick_spe = [[3000, 6000, 9000, 12000], 
        ['3', '6', '9', '12']]

XLABEL = r'Photometric redshift ($z_{\rm B}$)'
YLABEL = r'Counts per deg$^2$ per $0.1~z_{\rm B}$ [$\times 10^{3}$]'
DENSITY = False

plotting.HistPlotFunc(outpath,
                paras, wgs, COLORs, LABELs,
                nbins, XRANGE, YRANGE=None,
                XLABEL=XLABEL, YLABEL=YLABEL,
                DENSITY=DENSITY, HISTTYPE='step', STACKED=False,
                TITLE=None, xtick_min_label=True, ytick_min_label=True,
                xtick_spe=xtick_spe, ytick_spe=ytick_spe,
                vlines=None, vline_styles=None, vline_colors=None, vline_labels=None, vline_widths=None,
                hlines=None, hline_styles=None, hline_colors=None, hline_labels=None, hline_widths=None,
                xlog=False, ylog=False,
                loc_legend='best', 
                font_size=font_size, usetex=True,
                cumulative=False, 
                LINEs=LINEs, LINEWs=LINEWs)
