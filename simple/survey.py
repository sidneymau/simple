#!/usr/bin/env python
"""
"""
__author__ = "Sidney Mau"

import glob
import os

import numpy as np

import ugali.utils.healpix

#----------------------------------------------------------

class Survey():
    def __init__(self, iterable=(), **kwargs):
        self.__dict__.update(iterable, **kwargs)

        self.mag_1 = self.mag.format(self.band_1.upper())
        self.mag_2 = self.mag.format(self.band_2.upper())
        self.mag_err_1 = self.mag_err.format(self.band_1.upper())
        self.mag_err_2 = self.mag_err.format(self.band_2.upper())
        self.mag_dered_1 = self.mag_dered.format(self.band_1.upper())
        self.mag_dered_2 = self.mag_dered.format(self.band_2.upper())

        self.cols = [self.basis_1, self.basis_2,
                     self.mag_1, self.mag_2,
                     self.mag_err_1, self.mag_err_2,
                     self.mag_dered_1, self.mag_dered_2]

    def get_data(self, pixels):
        data_array = []
        for pixel in pixels:
            inlist = glob.glob('{}/*_{:05d}.fits'.format(self.datadir, pixel))
            for infile in inlist:
                if not os.path.exists(infile):
                    continue
                data_array.append(fits.read(infile, columns=self.cols))
        data = np.concatenate(data_array)
        return(data)

    def get_fracdet(self):
        if (self.fracdet is not None) and (self.fracdet.lower().strip() != 'none') and (self.fracdet != ''):
            print('Reading fracdet map {} ...').format(self.fracdet)
            fracdet = ugali.utils.healpix.read_map(self.fracdet)
        else:
            print('No fracdet map specified ...')
            fracdet = None
        return(fracdet)
