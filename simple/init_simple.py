#!/usr/bin/env python
"""
"""
__author__ = "Sidney Mau"

# Python libraries
import os
import yaml
import argparse

#----------------------------------------------------------

config = {'band_1' : 'G',
          'band_2' : 'R',
          'band_3' : 'I',
          'survey': {'fracdet' : '/Users/mcnanna/Research/y6/y6a1_griz_o.4096_t.32768_coverfoot_EQU.fits'},
          'catalog' : {'dirname'   : '/Users/mcnanna/Research/y6/far-out/datafiles/y6a1_skim',
                       'nside'     : 32,
                       'mag_max'   : 99.,
                       'basis_1'   : 'RA',
                       'basis_2'   : 'DEC',
                       'mag'       : 'SOF_PSF_MAG_CORRECTED_{}',
                       'mag_err'   : 'SOF_PSF_MAG_ERR_{}',
                       'reddening' : None,
                       'quality'   : 'SOF_PSF_MAG_CORRECTED_R < 24.5 && SOF_PSF_MAG_CORRECTED_I < 24.25',
                       'stars'     : 'EXT_SOF >= 0 && EXT_SOF <= 2',
                       'galaxies'  : 'EXT_SOF > 2',
                       'other'     : 'colorcolor'},
          'isochrone': {'name'        : 'Bressan2012',
                        'survey'      : 'des',
                        'age'         : 12,
                        'metallicity' : 0.0001},
          'output' : {'results_dir' : 'results_dir',
                      'log_dir'     : 'log_dir',
                      'plot_dir'    : 'plot_dir'},
          'grid' : {'delta_x'   : 0.003,
                    'smoothing' : 1.0, 
                    'bin_edge'  : 8.0,
                    'grid_dir'   : '/Users/mcnanna/Research/y6/far-out/datafiles/y6a1_skim/cutcut_arrays'},
          'moduli' : {'start' : 22.5,
                      'end'   : 26.5,
                      'step'  : 0.5},
          'batch' : {'max_jobs' : 20}
         }

def init_dirs(cfg):
    results_dir = os.path.join(os.getcwd(), cfg['output']['results_dir'])
    if not os.path.exists(results_dir):
        os.mkdir(results_dir)

    log_dir = os.path.join(os.getcwd(), cfg['output']['log_dir'])
    if not os.path.exists(log_dir):
        os.mkdir(log_dir)

    plot_dir = os.path.join(os.getcwd(), cfg['output']['plot_dir'])
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)

if __name__ == '__main__':
    cfgfile = 'config.yaml'
    with open(os.path.join(os.getcwd(), cfgfile), 'w') as outfile:
        yaml.dump(config, outfile, default_flow_style=False)

    with open(cfgfile, 'r') as ymlfile:
        cfg = yaml.load(ymlfile)

    init_dirs(cfg)
