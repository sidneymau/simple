#!/usr/bin/env python
"""
"""
__author__ = "Sidney Mau"

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
