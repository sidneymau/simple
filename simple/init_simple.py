#!/usr/bin/env python
"""
"""
__author__ = "Sidney Mau"

# Python libraries
import os
import yaml
import argparse

# Simple Libraries
import simple.survey

#----------------------------------------------------------

config = {'survey': {'name' : 'des',
                     'fracdet' : None,
                     'datadir': '/data/des40.b/data/y6a1/gold/1.1/healpix',
                     'nside'   : 32,
                     'mag_max' : 24.5,
                     'basis_1' : 'RA',
                     'basis_2' : 'DEC',
                     'band_1'  : 'G',
                     'band_2'  : 'R',
                     'mag'     : 'SOF_PSF_MAG_{}',
                     'mag_err' : 'SOF_PSF_MAG_ERR_{}',
                     'quality' : 'SOF_PSF_MAG_G < 24.5',
                     'stars'   : 'WAVG_SPREAD_MODEL_G < 0.003 + SPREADERR_MODEL_G'},
          'iso': {'name'   : 'Bressan2012',
                  'survey' : 'des'},
          'output' : {'results_dir' : 'results_dir',
                      'log_dir'     : 'log_dir',
                      'save_dir'    : 'plot_dir'},
          'jobs': 20
         }

def init_dirs(survey):
    results_dir = os.path.join(os.getcwd(), survey.output['results_dir'])
    if not os.path.exists(results_dir):
        os.mkdir(results_dir)

    log_dir = os.path.join(os.getcwd(), survey.output['log_dir'])
    if not os.path.exists(log_dir):
        os.mkdir(log_dir)

    save_dir = os.path.join(os.getcwd(), survey.output['save_dir'])
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)

if __name__ == '__main__':
    with open(os.path.join(os.getcwd(), 'config.yaml'), 'w') as outfile:
        yaml.dump(config, outfile, default_flow_style=False)

    #with open('config.yaml', 'r') as ymlfile:
    #    cfg = yaml.load(ymlfile)
    #    survey = simple.survey.Survey(cfg)

    #init_dirs(survey)

    ##import pdb;pdb.set_trace()
