#!/usr/bin/env python
"""
"""
__author__ = "Sidney Mau"

# Python libraries
#import os
#import yaml

#----------------------------------------------------------

class Survey():
    def __init__(self, iterable=(), **kwargs):
        self.__dict__.update(iterable, **kwargs)

        #self.survey = cfg['survey']
        #self.nside = cfg['nside']
        #self.mag_max = cfg['mag_max']

        #self.datadir = cfg['datadir']

        #self.basis_1 = cfg['basis_1']
        #self.basis_2 = cfg['basis_2']

        #self.band_1 = cfg['band_1']
        #self.band_2 = cfg['band_2']

        #self.mag = cfg['mag']
        #self.mag_err = cfg['mag_err']
        #self.mag_dered = cfg['mag_dered']

        #self.fracdet = cfg['fracdet']

        #self.mode = cfg['mode']
        #self.sim_population = cfg['sim_population']
        #self.sim_dir = cfg['sim_dir']

        self.mag_1 = self.mag.format(self.band_1.upper())
        self.mag_2 = self.mag.format(self.band_2.upper())
        self.mag_err_1 = self.mag_err.format(self.band_1.upper())
        self.mag_err_2 = self.mag_err.format(self.band_2.upper())
        self.mag_dered_1 = self.mag_dered.format(self.band_1.upper())
        self.mag_dered_2 = self.mag_dered.format(self.band_2.upper())


# The following should go in some initial setup
#with open('config.yaml', 'r') as ymlfile:
#    cfg = yaml.load(ymlfile)
#    survey = simple.survey.Survey(cfg)

#def init_dir():
#    results_dir = os.path.join(os.getcwd(), cfg['output']['results_dir'])
#    if not os.path.exists(results_dir):
#        os.mkdir(results_dir)
#
#    log_dir = os.path.join(os.getcwd(), cfg['output']['log_dir'])
#    if not os.path.exists(log_dir):
#        os.mkdir(log_dir)
#
#    save_dir = os.path.join(os.getcwd(), cfg['output']['save_dir'])
#    if not os.path.exists(save_dir):
#        os.mkdir(save_dir)
