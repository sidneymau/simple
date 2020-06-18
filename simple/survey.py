#!/usr/bin/env python
"""
"""
__author__ = "Sidney Mau"

import glob
import os
import importlib

import numpy as np
import healpy as hp
import fitsio as fits
import scipy.ndimage

import ugali.utils.healpix
import ugali.utils.projector
import ugali.isochrone

import simple.isochrone

#------------------------------------------------------------------------------

class Survey():
    """
    Class to handle survey-specific parameters.
    """
    def __init__(self, iterable=(), **kwargs):
        self.__dict__.update(iterable, **kwargs)

        self.mag_1 = self.catalog['mag'].format(self.band_1.upper())
        self.mag_2 = self.catalog['mag'].format(self.band_2.upper())
        self.mag_err_1 = self.catalog['mag_err'].format(self.band_1.upper())
        self.mag_err_2 = self.catalog['mag_err'].format(self.band_2.upper())
        self.cols = [self.catalog['basis_1'], self.catalog['basis_2'],
                     self.mag_1, self.mag_2,
                     self.mag_err_1, self.mag_err_2]
        if self.band_3 is not None:
            self.mag_3 = self.catalog['mag'].format(self.band_3.upper())
            self.mag_err_3 = self.catalog['mag_err'].format(self.band_3.upper())
            self.cols += [self.mag_3, self.mag_err_3]
        else:
            self.mag_3 = None
            self.mag_err_3 = None

        if self.catalog['reddening']:
            self.reddening_1 = self.catalog['reddening'].format(self.band_1.upper())
            self.reddening_2 = self.catalog['reddening'].format(self.band_2.upper())
            self.cols.append(self.reddening_1)
            self.cols.append(self.reddening_2)
            if self.band_3 is not None:
                self.reddening_3 = self.catalog['reddening'].format(self.band_3.upper())
                self.cols.append(self.reddening_3)
            else:
                self.reddening_3 = None

        # Set defaults for some optional configs
        grid_defaults = {'delta_x':0.01, 'smoothing':2.0, 'bin_edge':8.0, 'grid_dir':None}
        if not (hasattr(self, 'grid') and self.grid is not None):
            self.grid = grid_defaults
        else:
            for key in grid_defaults.keys():
                if key not in self.grid.keys() or self.grid[key] is None:
                    self.grid[key] = grid_defaults[key]

        moduli_defaults = {'start':16, 'end':self.catalog['mag_max'], 'step':0.5}
        if not (hasattr(self, 'moduli') and self.moduli is not None):
            self.moduli = moduli_defaults
        else:
            for key in moduli_defaults.keys():
                if key not in self.moduli.keys() or self.moduli[key] is None:
                    self.moduli[key] = moduli_defaults[key]

        self.load_fracdet

    @property
    def load_fracdet(self):
        """
        Load-in the fracdet map if it exists.
        """
        if self.survey['fracdet']:
            print('Reading fracdet map {} ...'.format(self.survey['fracdet']))
            fracdet = ugali.utils.healpix.read_map(self.survey['fracdet'])
        else:
            print('No fracdet map specified ...')
            fracdet = None
        #return(fracdet)
        self.fracdet = fracdet

    def get_neighbors(self, ra, dec):
        """
        Return center healpixel and 8 nearest neighbors for a given ra, dec pair.
        """
        pix_select = ugali.utils.healpix.angToPix(self.catalog.nside, ra, dec)
        pix_neighbors = np.concatenate([[pix_select], hp.get_all_neighbours(self.catalog['nside'], pix_select)])
        return(pix_neighbors)

    def get_data(self, pixels, type='stars'):
        """
        Load-in and return data for a list of healpixels as a numpy array.
        """
        if self.catalog['stars'] is None:
            print("No star selection given. Assuming all objects are stars")
            star_sel = None
        else:
            if type == 'stars':
                star_sel = self.catalog['stars']
            elif type == 'galaxies':
                star_sel = self.catalog['galaxies']
            elif type == 'all':
                star_sel = None
        sel = ' && '.join([s for s in (star_sel, self.catalog['quality']) if s is not None])

        data_array = []
        for pixel in pixels:
            inlist = glob.glob('{}/*_{:05d}.fits'.format(self.catalog['dirname'], pixel))
            for infile in inlist:
                if not os.path.exists(infile):
                    continue
                with fits.FITS(infile,vstorage='object') as f:
                    if len(sel) > 0:
                        w = f[1].where(sel)
                        d = f[1][self.cols][w]
                    else:
                        d = f[1][self.cols][:]
                    if self.catalog['reddening']:
                        d[self.mag_1] -= d[self.reddening_1]
                        d[self.mag_2] -= d[self.reddening_2]
                        if self.band_3 is not None:
                            d[self.mag_3] -= d[self.reddening_3]
                        d = d[[name for name in d.dtype.names if name not in [self.reddening_1, self.reddening_2, self.reddening_3]]]
                    data_array.append(d)
        data = np.concatenate(data_array)
        if self.catalog['other'] is not None:
            for module in self.catalog['other'].split('&&'):
                try: # Python 2
                    other = importlib.import_module(module.strip())
                except: # Python 3
                    spec = importlib.util.spec_from_file_location(module, os.getcwd()+'/{}.py'.format(module))
                    other = importlib.util.module_from_spec(spec)
                    spec.loader.exec_module(other)
                data = data[other.sel(self, data)]
        return(data)

    def get_isochrone(self, distance_modulus=None):
        iso_type = type(ugali.isochrone.factory(name=self.isochrone['name']))
        iso = simple.isochrone.get_isochrone(iso_type, survey=self.isochrone['survey'])
        iso.age = self.isochrone['age']
        iso.metallicity = self.isochrone['metallicity']
        if distance_modulus is not None:
            iso.distance_modulus = distance_modulus
        return(iso)

#------------------------------------------------------------------------------

class Region():
    """
    Class to handle regions.
    """
    def __init__(self, survey, ra, dec):
        self.survey = survey
        self.nside = self.survey.catalog['nside']
        self.fracdet = self.survey.fracdet

        self.ra = ra
        self.dec = dec
        self.proj = ugali.utils.projector.Projector(self.ra, self.dec)
        self.pix_center = ugali.utils.healpix.angToPix(self.nside, self.ra, self.dec)
        self.pix_neighbors = np.concatenate([[self.pix_center], hp.get_all_neighbours(self.nside, self.pix_center)])
        self.characteristic_density = None

    def get_data(self, type='stars'):
        return(self.survey.get_data(self.pix_neighbors, type))

    def load_data(self, type='stars'):
        self.data = self.get_data(type)

    #def get_stars(self, data):
    #    return(self.survey.get_stars(data))

    #def get_galaxies(self, data):
    #    return(self.survey.get_galaxies(data))

    def compute_characteristic_density(self, data):
        """
        Compute the characteristic density of a region
        Convlve the field and find overdensity peaks
        """

        cut_magnitude_threshold = (data[self.survey.mag_1] < self.survey.catalog['mag_max'])
    
        x, y = self.proj.sphereToImage(data[self.survey.catalog['basis_1']][cut_magnitude_threshold], data[self.survey.catalog['basis_2']][cut_magnitude_threshold]) # Trimmed magnitude range for hotspot finding
    
        delta_x_coverage = 0.1
        area_coverage = (delta_x_coverage)**2
        bins_coverage = np.arange(-5., 5. + 1.e-10, delta_x_coverage)
        h_coverage = np.histogram2d(x, y, bins=[bins_coverage, bins_coverage])[0]
        #h_goodcoverage = np.histogram2d(x[cut_goodcoverage], y[cut_goodcoverage], bins=[bins_coverage, bins_coverage])[0]
        h_goodcoverage = np.histogram2d(x, y, bins=[bins_coverage, bins_coverage])[0]
    
        n_goodcoverage = h_coverage[h_goodcoverage > 0].flatten()
    
        #characteristic_density = np.mean(n_goodcoverage) / area_coverage # per square degree
        characteristic_density = np.median(n_goodcoverage) / area_coverage # per square degree
        print('Characteristic density = {:0.1f} deg^-2'.format(characteristic_density))
    
        # Use pixels with fracdet ~1.0 to estimate the characteristic density
        if self.fracdet is not None:
            fracdet_zero = np.tile(0., len(self.fracdet))
            cut = (self.fracdet != hp.UNSEEN)
            fracdet_zero[cut] = self.fracdet[cut]
    
            nside_fracdet = hp.npix2nside(len(self.fracdet))
            
            subpix_region_array = []
            for pix in np.unique(ugali.utils.healpix.angToPix(self.nside, data[self.survey.catalog['basis_1']], data[self.survey.catalog['basis_2']])):
                subpix_region_array.append(ugali.utils.healpix.subpixel(self.pix_center, self.nside, nside_fracdet))
            subpix_region_array = np.concatenate(subpix_region_array)
    
            # Compute mean fracdet in the region so that this is available as a correction factor
            cut = (self.fracdet[subpix_region_array] != hp.UNSEEN)
            mean_fracdet = np.mean(self.fracdet[subpix_region_array[cut]])
    
            # Correct the characteristic density by the mean fracdet value
            characteristic_density /= mean_fracdet 
            print('Characteristic density (fracdet corrected) = {:0.1f} deg^-2'.format(characteristic_density))
    
        return(characteristic_density)

    def set_characteristic_density(self, data):
        d = self.compute_characteristic_density(data)
        self.characteristic_density = d
    
    def characteristic_density_local(self, data, x_peak, y_peak, angsep_peak):
        """
        Compute the local characteristic density of a region
        """
    
        cut_magnitude_threshold = (data[self.survey.mag_1] < self.survey.catalog['mag_max'])

        if self.characteristic_density is None:
            characteristic_density = self.compute_characteristic_density(data)
        else:
            characteristic_density = self.characteristic_density
    
        x, y = self.proj.sphereToImage(data[self.survey.catalog['basis_1']][cut_magnitude_threshold], data[self.survey.catalog['basis_2']][cut_magnitude_threshold]) # Trimmed magnitude range for hotspot finding
        #x_full, y_full = proj.sphereToImage(data[basis_1], data[basis_2]) # If we want to use full magnitude range for significance evaluation
        inner, outer = 0.3, 0.5 
    
        # If fracdet map is available, use that information to either compute local density,
        # or in regions of spotty coverage, use the typical density of the region
        if self.fracdet is not None:
            # The following is copied from how it's used in compute_char_density
            fracdet_zero = np.tile(0., len(self.fracdet))
            cut = (self.fracdet != hp.UNSEEN)
            fracdet_zero[cut] = self.fracdet[cut]
    
            nside_fracdet = hp.npix2nside(len(self.fracdet))
            
            subpix_region_array = []
            for pix in np.unique(ugali.utils.healpix.angToPix(self.nside, data[self.survey.catalog['basis_1']], data[self.survey.catalog['basis_2']])):
                subpix_region_array.append(ugali.utils.healpix.subpixel(self.pix_center, self.nside, nside_fracdet))
            subpix_region_array = np.concatenate(subpix_region_array)
    
            # Compute mean fracdet in the region so that this is available as a correction factor
            cut = (self.fracdet[subpix_region_array] != hp.UNSEEN)
            mean_fracdet = np.mean(self.fracdet[subpix_region_array[cut]])
    
            subpix_region_array = subpix_region_array[self.fracdet[subpix_region_array] > 0.99]
            subpix = ugali.utils.healpix.angToPix(nside_fracdet, 
                                                  data[self.survey.catalog['basis_1']][cut_magnitude_threshold], 
                                                  data[self.survey.catalog['basis_2']][cut_magnitude_threshold])
    
            # This is where the local computation begins
            ra_peak, dec_peak = self.proj.imageToSphere(x_peak, y_peak)
            subpix_all = ugali.utils.healpix.angToDisc(nside_fracdet, ra_peak, dec_peak, outer)
            subpix_inner = ugali.utils.healpix.angToDisc(nside_fracdet, ra_peak, dec_peak, inner)
            subpix_annulus = subpix_all[~np.in1d(subpix_all, subpix_inner)]
            mean_fracdet = np.mean(fracdet_zero[subpix_annulus])
            print('mean_fracdet {}'.format(mean_fracdet))
            if mean_fracdet < 0.5:
                characteristic_density_local = characteristic_density
                print('Characteristic density local (baseline) = {:0.1f} deg^-2 = {:0.3f} arcmin^-2'.format(characteristic_density_local, characteristic_density_local / 60.**2))
            else:
                # Check pixels in annulus with complete coverage
                subpix_annulus_region = np.intersect1d(subpix_region_array, subpix_annulus)
                print('{} percent pixels with complete coverage'.format(float(len(subpix_annulus_region)) / len(subpix_annulus)))
                if (float(len(subpix_annulus_region)) / len(subpix_annulus)) < 0.25:
                    characteristic_density_local = characteristic_density
                    print('Characteristic density local (spotty)  = {:0.1f} deg^-2 = {:0.3f} arcmin^-2'.format(characteristic_density_local, characteristic_density_local / 60.**2))
                else:
                    characteristic_density_local = float(np.sum(np.in1d(subpix, subpix_annulus_region))) \
                                                   / (hp.nside2pixarea(nside_fracdet, degrees=True) * len(subpix_annulus_region)) # deg^-2
                    print('Characteristic density local (cleaned up) = {:0.1f} deg^-2 = {:0.3f} arcmin^-2'.format(characteristic_density_local, characteristic_density_local / 60.**2))
        else:
            # If not good azimuthal coverage, revert to broad region density
            cut_annulus = (angsep_peak > inner) & (angsep_peak < outer)
            #phi = np.degrees(np.arctan2(y_full[cut_annulus] - y_peak, x_full[cut_annulus] - x_peak)) # Use full magnitude range, NOT TESTED!!!
            phi = np.degrees(np.arctan2(y[cut_annulus] - y_peak, x[cut_annulus] - x_peak)) # Impose magnitude threshold
            h = np.histogram(phi, bins=np.linspace(-180., 180., 13))[0]
            if np.sum(h > 0) < 10 or np.sum(h > 0.5 * np.median(h)) < 10:
                characteristic_density_local = characteristic_density
            else:
                # Compute the local characteristic density
                area_field = np.pi * (outer**2 - inner**2)
                n_field = np.sum((angsep_peak > inner) & (angsep_peak < outer))
                characteristic_density_local = n_field / area_field
            print('Characteristic density local = {:0.1f} deg^-2 = {:0.3f} arcmin^-2'.format(characteristic_density_local, characteristic_density_local / 60.**2))
    
        return(characteristic_density_local)

    def find_peaks(self, data, distance_modulus):
        """
        Convolve field to find characteristic density and peaks within the selected pixel
        """

        # convolve field and find peaks
        cut_magnitude_threshold = (data[self.survey.mag_1] < self.survey.catalog['mag_max'])

        if self.characteristic_density is None:
            characteristic_density = self.compute_characteristic_density(data)
        else:
            characteristic_density = self.characteristic_density
    
        x, y = self.proj.sphereToImage(data[self.survey.catalog['basis_1']][cut_magnitude_threshold], data[self.survey.catalog['basis_2']][cut_magnitude_threshold]) # Trimmed magnitude range for hotspot finding
        #x_full, y_full = proj.sphereToImage(data[basis_1], data[basis_2]) # If we want to use full magnitude range for significance evaluation
        area = self.survey.grid['delta_x']**2
        bins = np.arange(-self.survey.grid['bin_edge'], self.survey.grid['bin_edge'] + 1.e-10, self.survey.grid['delta_x'])
        centers = 0.5 * (bins[0: -1] + bins[1:])
        yy, xx = np.meshgrid(centers, centers)
    
        h = np.histogram2d(x, y, bins=[bins, bins])[0]
        h_g = scipy.ndimage.filters.gaussian_filter(h, (self.survey.grid['smoothing']/60.) / self.survey.grid['delta_x'])
        
        if self.survey.grid['grid_dir'] is None:
            rara, decdec = self.proj.imageToSphere(xx.flatten(), yy.flatten())
            cutcut = (ugali.utils.healpix.angToPix(self.nside, rara, decdec) == self.pix_center).reshape(xx.shape)
        else:
            cutcut_file = glob.glob('{}/*_{}.npz'.format(self.survey.grid['grid_dir'], self.pix_center))[0]
            cutcut = np.load(cutcut_file)['arr_0']

        # Course first to find cutoff region
        for factor in [5., 4., 3., 2., 1.]:
            threshold_density = area * characteristic_density * factor
            h_region, n_region = scipy.ndimage.measurements.label((h_g * cutcut) > (threshold_density))
            #print 'factor', factor, n_region, n_region < 10
            if n_region >= 10:
                break
        # Fine to find exact threshold density
        factor_array = np.arange(factor+0.05, 5., 0.05)
        for factor in factor_array:
            threshold_density = area * characteristic_density * factor
            h_region, n_region = scipy.ndimage.measurements.label((h_g * cutcut) > (threshold_density))
            #print 'factor', factor, n_region, n_region < 10
            if n_region < 10:
                break
    
        h_region = np.ma.array(h_region, mask=(h_region < 1))
    
        x_peak_array = []
        y_peak_array = []
        angsep_peak_array = []
    
        for index in range(1, n_region + 1): # loop over peaks
            index_peak = np.argmax(h_g * (h_region == index))
            x_peak, y_peak = xx.flatten()[index_peak], yy.flatten()[index_peak]
            #print index, np.max(h_g * (h_region == index))
            
            #angsep_peak = np.sqrt((x_full - x_peak)**2 + (y_full - y_peak)**2) # Use full magnitude range, NOT TESTED!!!
            angsep_peak = np.sqrt((x-x_peak)**2 + (y-y_peak)**2) # Impose magnitude threshold
    
            x_peak_array.append(x_peak)
            y_peak_array.append(y_peak)
            angsep_peak_array.append(angsep_peak)
        
        return x_peak_array, y_peak_array, angsep_peak_array
    
    def fit_aperture(self, data, distance_modulus, x_peak, y_peak, angsep_peak):
        """
        Fit aperture by varing radius and computing the significance
        """

        characteristic_density_local = self.characteristic_density_local(data, x_peak, y_peak, angsep_peak)
    
        ra_peak_array = []
        dec_peak_array = []
        r_peak_array = []
        sig_peak_array = []
        distance_modulus_array = []
        n_obs_peak_array = []
        n_obs_half_peak_array = []
        n_model_peak_array = []
    
        size_array = np.arange(0.01, 0.3, 0.01)
        #size_array = np.concatenate((np.arange(0.003, 0.01, 0.001), np.arange(0.01, 0.3, 0.01)))
        sig_array = np.tile(0., len(size_array))
        
        size_array_zero = np.concatenate([[0.], size_array])
        area_array = np.pi * (size_array_zero[1:]**2 - size_array_zero[0:-1]**2)
    
        n_obs_array = np.tile(0, len(size_array))
        n_model_array = np.tile(0., len(size_array))
        for ii in range(0, len(size_array)):
            n_obs = np.sum(angsep_peak < size_array[ii])
            n_model = characteristic_density_local * (np.pi * size_array[ii]**2)
            sig_array[ii] = np.clip(scipy.stats.norm.isf(scipy.stats.poisson.sf(n_obs, n_model)), 0., 37.5) # Clip at 37.5
            n_obs_array[ii] = n_obs
            n_model_array[ii] = n_model
    
        ra_peak, dec_peak = self.proj.imageToSphere(x_peak, y_peak)
    
        index_peak = np.argmax(sig_array)
        r_peak = size_array[index_peak]
        #if np.max(sig_array) >= 37.5:
        #    r_peak = 0.5
        n_obs_peak = n_obs_array[index_peak]
        n_model_peak = n_model_array[index_peak]
        n_obs_half_peak = np.sum(angsep_peak < (0.5 * r_peak))
    
        # Compile resilts
        print('Candidate: x_peak: {:12.3f}, y_peak: {:12.3f}, r_peak: {:12.3f}, sig: {:12.3f}, ra_peak: {:12.3f}, dec_peak: {:12.3f}'.format(x_peak, y_peak, r_peak, np.max(sig_array), ra_peak, dec_peak))
        ra_peak_array.append(ra_peak)
        dec_peak_array.append(dec_peak)
        r_peak_array.append(r_peak)
        #sig_peak_array.append(np.max(sig_array))
        sig_peak_array.append(sig_array[index_peak])
        distance_modulus_array.append(distance_modulus)
        n_obs_peak_array.append(n_obs_peak)
        n_obs_half_peak_array.append(n_obs_half_peak)
        n_model_peak_array.append(n_model_peak)
    
        return ra_peak_array, dec_peak_array, r_peak_array, sig_peak_array, distance_modulus_array, n_obs_peak_array, n_obs_half_peak_array, n_model_peak_array
