#!/usr/bin/env python
"""
"""
__author__ = "Sidney Mau"

import sys
import os
import glob
import yaml
import argparse
import numpy as np
import healpy as hp
import fitsio as fits

import ugali.utils.healpix
import ugali.utils.projector

import simple.filters
import simple.simple_utils
import simple.survey

#------------------------------------------------------------------------------

def search_by_distance(survey, region, data, distance_modulus):
    """
    Idea: 
    Send a data extension that goes to faint magnitudes, e.g., g < 24.
    Use the whole region to identify hotspots using a slightly brighter 
    magnitude threshold, e.g., g < 23, so not susceptible to variations 
    in depth. Then compute the local field density using a small annulus 
    around each individual hotspot, e.g., radius 0.3 to 0.5 deg.
    """

    print('Distance = {:0.1f} kpc (m-M = {:0.1f})').format(ugali.utils.projector.distanceModulusToDistance(distance_modulus), distance_modulus)

    iso = ugali.isochrone.factory(name=survey.isoname, survey=survey.isosurvey, band_1=survey.band_1.lower(), band_2=survey.band_2.lower())
    iso.age = 12.
    iso.metallicity = 0.0001
    iso.distance_modulus = distance_modulus

    cut = simple_utils.cut_isochrone_path(data[survey.mag_dered_1], data[survey.mag_dered_2], data[survey.mag_err_1], data[survey.mag_err_2], iso, radius=0.1)
    data = data[cut]

    print('{} objects left after isochrone cut...').format(len(data))

    if (len(data) == 0):
        return [], [], [], [], [], [], [], []

    characteristic_density = region.characteristic_density(data)

    ra_peak_array = []
    dec_peak_array = []
    r_peak_array = []
    sig_peak_array = []
    distance_modulus_array = []
    n_obs_peak_array = []
    n_obs_half_peak_array = []
    n_model_peak_array = []

    x_peak_array, y_peak_array, angsep_peak_array = region.find_peaks(data, distance_modulus)

    for x_peak, y_peak, angsep_peak in itertools.izip(x_peak_array, y_peak_array, angsep_peak_array):
        characteristic_density_local = region.characteristic_density_local(data, x_peak, y_peak, angsep_peak)
        # Aperture fitting
        print('Fitting aperture to hotspot...')
        ra_peaks, dec_peaks, r_peaks, sig_peaks, distance_moduli, n_obs_peaks, n_obs_half_peaks, n_model_peaks = region.fit_aperture(distance_modulus, x_peak, y_peak, angsep_peak)
        
        ra_peak_array.append(ra_peaks)
        dec_peak_array.append(dec_peaks)
        r_peak_array.append(r_peaks)
        sig_peak_array.append(sig_peaks)
        distance_modulus_array.append(distance_moduli)
        n_obs_peak_array.append(n_obs_peaks)
        n_obs_half_peak_array.append(n_obs_half_peaks)
        n_model_peak_array.append(n_model_peaks)

    ra_peak_array = np.concatenate(ra_peak_array)
    dec_peak_array = np.concatenate(dec_peak_array)
    r_peak_array = np.concatenate(r_peak_array)
    sig_peak_array = np.concatenate(sig_peak_array)
    distance_modulus_array = np.concatenate(distance_modulus_array)
    n_obs_peak_array = np.concatenate(n_obs_peak_array)
    n_obs_half_peak_array = np.concatenate(n_obs_half_peak_array)
    n_model_peak_array = np.concatenate(n_model_peak_array)

    return ra_peak_array, dec_peak_array, r_peak_array, sig_peak_array, distance_modulus_array, n_obs_peak_array, n_obs_half_peak_array, n_model_peak_array

#------------------------------------------------------------------------------

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

    survey.load_fracdet

    #--------------------------------------------------------------------------

    region = simple.survey.Region(survey, args['ra'], args['dec'])
    print('Search coordinates: (RA, Dec) = ({:0.2f}, {:0.2f})').format(region.ra, region.dec)
    print('Search healpixel: {} (nside = {})').format(region.pix_center, region.nside)
    print('Healpixels: {}'.format(region.pix_neighbors))

    #--------------------------------------------------------------------------

    data = region.get_data()
    # assume data is cut for quality and dereddened during skim
    print('Found {} objects').format(len(data))
    if (len(data) == 0):
        print('Ending search.')
        exit()

    star_sel = region.get_stars(data)
    galaxy_sel = region.get_galaxies(data)
    print('Found {} stars').format(sum(star_sel))
    print('Found {} galaxies').format(sum(galaxy_sel))
    if (sum(star_sel) == 0):
        print('Ending search.')
        exit()

    stars = data[star_sel]

    #--------------------------------------------------------------------------

    distance_modulus_search_array = np.arange(16., survey.mag_max, 0.5)

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
        ra_peaks, dec_peaks, r_peaks, sig_peaks, dist_moduli, n_obs_peaks, n_obs_half_peaks, n_model_peaks = simple.simple_utils.search_by_distance(survey, region, data[stars], distance_modulus)
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
        simple.simple_utils.write_output(results_dir, survey.nside, pix_nside_select, ra_peak_array, dec_peak_array,
                                         r_peak_array, distance_modulus_array, 
                                         n_obs_peak_array, n_obs_half_peak_array, n_model_peak_array, 
                                         sig_peak_array, mc_source_id_array, mode, outfile)
    else:
        print('No significant hotspots found.')
        nan_array = [np.nan]
        simple.simple_utils.write_output(results_dir, survey.nside, pix_nside_select,
                                         nan_array, nan_array, nan_array, nan_array, 
                                         nan_array, nan_array, nan_array, nan_array,
                                         [mc_source_id], mode, outfile)
