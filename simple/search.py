#!/usr/bin/env python
"""
"""
__author__ = "Sidney Mau"

# Python libraries
import sys
import os
import glob
import yaml
from matplotlib import mlab
import numpy as np
import healpy as hp
import astropy.io.fits as pyfits
import fitsio as fits
import argparse

# Ugali libraries
import ugali.utils.healpix
import ugali.utils.projector

# Simple binner modules
import simple.filters
import simple.simple_utils
import simple.survey

#----------------------------------------------------------

if __name__ == '__main__':
    # Construct argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--infile',type=str,required=False,
                        help='Input file')
    parser.add_argument('--outfile',type=str,required=False,
                        help='Output file')
    parser.add_argument('--ra',type=float,required=False,
                        help='Target RA')
    parser.add_argument('--dec',type=float,required=False,
                        help='Target DEC')
    args = vars(parser.parse_args())

    with open('defaults.yaml', 'r') as ymlfile:
        cfg = yaml.load(ymlfile)
        survey = simple.survey.Survey(cfg)

    #-------------------------------

    ra = args['ra']
    dec = args['dec']
    print('Search coordinates: (RA, Dec) = ({:0.2f}, {:0.2f})').format(ra, dec)

    neighbors = survey.get_neighbors(ra, dec)
    print('Healpixels: {}'.format(neighbors))

    #-------------------------------

    data = survey.get_data(neighbors)
    # assume data is cut for quality and dereddened during skim
    print('Found {} objects').format(len(data))
    if (len(data) == 0):
        print('Ending search.')
        exit()

    star_sel = survey.get_stars(data)
    galaxy_sel = survey.get_galaxies(data)
    print('Found {} stars').format(sum(star_sel))
    print('Found {} galaxies').format(sum(galaxy_sel))
    if (sum(star_sel) == 0):
        print('Ending search.')
        exit()

    stars = data[star_sel]

    #-------------------------------

    fracdet = survey.get_fracdet()

    #-------------------------------

    distance_modulus_search_array = np.arange(16., survey.mag_max, 0.5)

############################################################

ra_peak_array = []
dec_peak_array = [] 
r_peak_array = []
sig_peak_array = []
distance_modulus_array = []
mc_source_id_array = []
n_obs_peak_array = []
n_obs_half_peak_array = []
n_model_peak_array = []

for distance_modulus in distance_modulus_search_array:
    ra_peaks, dec_peaks, r_peaks, sig_peaks, dist_moduli, n_obs_peaks, n_obs_half_peaks, n_model_peaks = simple.simple_utils.search_by_distance(nside, data, distance_modulus, pix_nside_select, ra_select, dec_select, mag_max, fracdet)
    ra_peak_array.append(ra_peaks)
    dec_peak_array.append(dec_peaks)
    r_peak_array.append(r_peaks)
    sig_peak_array.append(sig_peaks)
    distance_modulus_array.append(dist_moduli)
    n_obs_peak_array.append(n_obs_peaks)
    n_obs_half_peak_array.append(n_obs_half_peaks)
    n_model_peak_array.append(n_model_peaks)
    mc_source_id_array.append(np.tile(0, len(sig_peaks)))

ra_peak_array = np.concatenate(ra_peak_array)
dec_peak_array = np.concatenate(dec_peak_array)
r_peak_array = np.concatenate(r_peak_array)
sig_peak_array = np.concatenate(sig_peak_array)
distance_modulus_array = np.concatenate(distance_modulus_array)
n_obs_peak_array = np.concatenate(n_obs_peak_array)
n_obs_half_peak_array = np.concatenate(n_obs_half_peak_array)
n_model_peak_array = np.concatenate(n_model_peak_array)
mc_source_id_array = np.concatenate(mc_source_id_array)

# Sort peaks according to significance
index_sort = np.argsort(sig_peak_array)[::-1]
ra_peak_array = ra_peak_array[index_sort]
dec_peak_array = dec_peak_array[index_sort]
r_peak_array = r_peak_array[index_sort]
sig_peak_array = sig_peak_array[index_sort]
distance_modulus_array = distance_modulus_array[index_sort]
n_obs_peak_array = n_obs_peak_array[index_sort]
n_obs_half_peak_array = n_obs_half_peak_array[index_sort]
n_model_peak_array = n_model_peak_array[index_sort]
mc_source_id_array = mc_source_id_array[index_sort]

# Collect overlapping peaks
for ii in range(0, len(sig_peak_array)):
    if sig_peak_array[ii] < 0:
        continue
    angsep = ugali.utils.projector.angsep(ra_peak_array[ii], dec_peak_array[ii], ra_peak_array, dec_peak_array)
    sig_peak_array[(angsep < r_peak_array[ii]) & (np.arange(len(sig_peak_array)) > ii)] = -1.
    #sig_peak_array[(angsep < 0.5) & (np.arange(len(sig_peak_array)) > ii)] = -1. # 0.5 deg radius

if (mode == 0):
    # Prune the list of peaks
    ra_peak_array = ra_peak_array[sig_peak_array > 0.]
    dec_peak_array = dec_peak_array[sig_peak_array > 0.]
    r_peak_array = r_peak_array[sig_peak_array > 0.]
    distance_modulus_array = distance_modulus_array[sig_peak_array > 0.]
    n_obs_peak_array = n_obs_peak_array[sig_peak_array > 0.]
    n_obs_half_peak_array = n_obs_half_peak_array[sig_peak_array > 0.]
    n_model_peak_array = n_model_peak_array[sig_peak_array > 0.]
    mc_source_id_array = mc_source_id_array[sig_peak_array > 0.]
    sig_peak_array = sig_peak_array[sig_peak_array > 0.] # Update the sig_peak_array last!

for ii in range(0, len(sig_peak_array)):
    print('{:0.2f} sigma; (RA, Dec, d) = ({:0.2f}, {:0.2f}); r = {:0.2f} deg; d = {:0.1f}, mu = {:0.2f} mag), mc_source_id: {:0.2f}'.format(sig_peak_array[ii], 
                 ra_peak_array[ii], 
                 dec_peak_array[ii], 
                 r_peak_array[ii],
                 ugali.utils.projector.distanceModulusToDistance(distance_modulus_array[ii]),
                 distance_modulus_array[ii],
                 mc_source_id_array[ii]))

# Write output
if (len(sig_peak_array) > 0):
    simple.simple_utils.write_output(results_dir, nside, pix_nside_select, ra_peak_array, dec_peak_array, r_peak_array, distance_modulus_array, 
                             n_obs_peak_array, n_obs_half_peak_array, n_model_peak_array, 
                             sig_peak_array, mc_source_id_array, mode, outfile)
else:
    print('No significant hotspots found.')
    nan_array = [np.nan]
    simple.simple_utils.write_output(results_dir, nside, pix_nside_select,
                             nan_array, nan_array, nan_array, nan_array, 
                             nan_array, nan_array, nan_array, nan_array,
                             [mc_source_id], mode, outfile)
