#!/usr/bin/env python
"""
"""
__author__ = "Sidney Mau"

# Python libraries
import os
import yaml
import argparse

#----------------------------------------------------------

config = {'survey': {'name'    : 'des',
                     'fracdet' : None},
          'band_1' : 'G',
          'band_2' : 'R',
          'band_3' : 'I',
          'catalog' : {'dirname'   : '/Users/mcnanna/Research/y6/far-out/datafiles/y6a1_skim',
                       'nside'     : 32,
                       'mag_max'   : 24.5,
                       'basis_1'   : 'RA',
                       'basis_2'   : 'DEC',
                       'mag'       : 'SOF_PSF_MAG_CORRECTED_{}',
                       'mag_err'   : 'SOF_PSF_MAG_ERR_{}',
                       'reddening' : None,
                       'quality'   : 'SOF_PSF_MAG_CORRECTED_R < 24.5 && SOF_PSF_MAG_CORRECTED_I < 24.25',
                       'stars'     : 'EXT_SOF >= 0 && EXT_SOF <= 2',
                       'other'     : 'colorcolor'},
          'isochrone': {'name'        : 'Bressan2012',
                        'survey'      : 'des',
                        'age'         : 12,
                        'metallicity' : 0.0001},
          'output' : {'results_dir' : 'results_dir',
                      'log_dir'     : 'log_dir',
                      'save_dir'    : 'plot_dir'},
          'batch' : {'max_jobs' : 20}
         }

def init_dirs(cfg):
    results_dir = os.path.join(os.getcwd(), cfg['output']['results_dir'])
    if not os.path.exists(results_dir):
        os.mkdir(results_dir)

    log_dir = os.path.join(os.getcwd(), cfg['output']['log_dir'])
    if not os.path.exists(log_dir):
        os.mkdir(log_dir)

    save_dir = os.path.join(os.getcwd(), cfg['output']['save_dir'])
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)

if __name__ == '__main__':
    cfgfile = 'config.yaml'
    with open(os.path.join(os.getcwd(), cfgfile), 'w') as outfile:
        yaml.dump(config, outfile, default_flow_style=False)

    with open(cfgfile, 'r') as ymlfile:
        cfg = yaml.load(ymlfile)

    init_dirs(cfg)
