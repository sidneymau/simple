#!/usr/bin/env python
"""
"""
__author__ = "Sidney Mau"

import glob
import os

import numpy as np
import healpy as hp

import ugali.utils.healpix

#------------------------------------------------------------------------------

class Survey():
    """
    Class to handle survey-specific parameters.
    """
    def __init__(self, iterable=(), **kwargs):
        self.__dict__.update(iterable, **kwargs)

        self.mag_1 = self.mag.format(self.band_1.upper())
        self.mag_2 = self.mag.format(self.band_2.upper())
        self.mag_err_1 = self.mag_err.format(self.band_1.upper())
        self.mag_err_2 = self.mag_err.format(self.band_2.upper())
        self.mag_dered_1 = self.mag_dered.format(self.band_1.upper())
        self.mag_dered_2 = self.mag_dered.format(self.band_2.upper())

        self.spread_model = 'WAVG_SPREAD_MODEL_{}'.format(self.band_1.upper())
        self.spreaderr_model = 'SPREADERR_MODEL_{}'.format(self.band_1.upper())

        self.cols = [self.basis_1, self.basis_2,
                     self.mag_1, self.mag_2,
                     self.mag_err_1, self.mag_err_2,
                     self.mag_dered_1, self.mag_dered_2,
                     self.spread_model, self.spreaderr_model]

        self.fracdet = None

    @property
    def load_fracdet(self):
        """
        Load-in the fracdet map if it exists.
        """
        if (self.fracdet is not None) and (self.fracdet.lower().strip() != 'none') and (self.fracdet != ''):
            print('Reading fracdet map {} ...').format(self.fracdet)
            fracdet = ugali.utils.healpix.read_map(self.fracdet)
        else:
            print('No fracdet map specified ...')
            fracdet = None
        #return(fracdet)
        self.fracdet = fracdet

    def get_neighbors(self, ra, dec):
        """
        Return center healpixel and 8 nearest neighbors for a given ra, dec pair.
        """
        pix_select = ugali.utils.healpix.angToPix(self.nside, ra, dec)
        pix_neighbors = np.concatenate([[pix_select], hp.get_all_neighbours(self.nside, pix_select)])
        return(pix_neighbors)

    def get_data(self, pixels):
        """
        Load-in and return data for a list of healpixels as a numpy array.
        """
        data_array = []
        for pixel in pixels:
            inlist = glob.glob('{}/*_{:05d}.fits'.format(self.datadir, pixel))
            for infile in inlist:
                if not os.path.exists(infile):
                    continue
                data_array.append(fits.read(infile, columns=self.cols))
        data = np.concatenate(data_array)
        return(data)

    def get_stars(self, data):
        """
        Return boolean list of star-like objects.
        """
        sel = (data['WAVG_SPREAD_MODEL_G'] < 0.003 + data['SPREADERR_MODEL_G'])
        return(sel)

    def get_galaxies(self, data):
        """
        Return boolean list of galaxy-like objects.
        """
        sel = (data['WAVG_SPREAD_MODEL_G'] > 0.003 + data['SPREADERR_MODEL_G'])
        return(sel)

class Region():
    """
    Class to handle regions.
    """
    def __init__(self, survey, ra, dec):
        self.survey = survey
        self.nside = self.survey.nside
        self.fracdet = self.survey.fracdet

        self.ra = ra
        self.dec = dec
        self.pix_center = ugali.utils.healpix.angToPix(self.nside, self.ra, self.dec)
        self.pix_neighbors = np.concatenate([[self.pix_center], hp.get_all_neighbours(self.nside, self.pix_center)])

    def get_data(self):
        return(self.survey.get_data(self.pix_neighbors))

    def get_stars(self, data):
        return(self.survey.get_stars(data))

    def get_galaxies(self, data):
        return(self.survey.get_galaxies(data))
