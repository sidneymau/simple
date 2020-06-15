#i!/usr/bin/env python
"""
"""
__author__ = "Sidney Mau"

import os
import yaml
import argparse
import itertools
import numpy as np
import scipy.interpolate

import ugali.utils.projector

import simple.survey

#------------------------------------------------------------------------------

def write_output(results_dir, nside, pix_nside_select, ra_peak_array, dec_peak_array, r_peak_array, distance_modulus_array, 
                n_obs_peak_array, n_obs_half_peak_array, n_model_peak_array, 
                sig_peak_array, mc_source_id_array, mode, outfile):
    #writer = open(outfile, 'a') # append if exists
    #for ii in range(0, len(sig_peak_array)):
    #    # SIG, RA, DEC, MODULUS, r, n_obs, n_model, mc_source_id
    #    writer.write('{:10.3f}, {:10.3f}, {:10.3f}, {:10.3f}, {:10.3f}, {:10.3f}, {:10.3f}, {:10.3f}, {:10.3f}\n'.format(sig_peak_array[ii], 
    #                                                                                                                     ra_peak_array[ii], 
    #                                                                                                                     dec_peak_array[ii], 
    #                                                                                                                     distance_modulus_array[ii], 
    #                                                                                                                     r_peak_array[ii],
    #                                                                                                                     n_obs_peak_array[ii],
    #                                                                                                                     n_obs_half_peak_array[ii],
    #                                                                                                                     n_model_peak_array[ii],
    #                                                                                                                     mc_source_id_array[ii]))
    #f = open(outfile, 'ab')
    #np.savetxt(f, arr, delimiter=',')
    #f.close()
    data = [tuple(row) for row in np.stack([sig_peak_array, ra_peak_array, dec_peak_array, distance_modulus_array, r_peak_array, n_obs_peak_array, n_obs_half_peak_array, n_model_peak_array, mc_source_id_array], axis=-1)]
    arr = np.array(data, dtype=[('SIG', float), ('RA', float), ('DEC', float), ('MODULUS', float), ('R', float), ('N_OBS', float), ('N_OBS_HALF', float), ('N_MODEL', float), ('MC_SOURCE_ID', int)])
    np.save('{}/{}'.format(results_dir, outfile), arr)


def search_by_distance(survey, region, distance_modulus):
    """
    Idea: 
    Send a data extension that goes to faint magnitudes, e.g., g < 24.
    Use the whole region to identify hotspots using a slightly brighter 
    magnitude threshold, e.g., g < 23, so not susceptible to variations 
    in depth. Then compute the local field density using a small annulus 
    around each individual hotspot, e.g., radius 0.3 to 0.5 deg.
    """

    print('Distance = {:0.1f} kpc (m-M = {:0.1f})'.format(ugali.utils.projector.distanceModulusToDistance(distance_modulus), distance_modulus))

    iso = survey.get_isochrone(distance_modulus)
    cut = iso.cut_separation(survey.band_1.lower(), survey.band_2.lower(), region.data[survey.mag_1], region.data[survey.mag_2], region.data[survey.mag_err_1], region.data[survey.mag_err_2], radius=0.1)
    if survey.band_3 is not None:
        cut &= iso.cut_separation(survey.band_2.lower(), survey.band_3.lower(), region.data[survey.mag_2], region.data[survey.mag_3], region.data[survey.mag_err_2], region.data[survey.mag_err_3], radius=0.1)
    data = region.data[cut]

    print('{} objects left after isochrone cut...'.format(len(data)))

    if (len(data) == 0):
        return [], [], [], [], [], [], [], []

    region.set_characteristic_density(data)

    ra_peak_array = []
    dec_peak_array = []
    r_peak_array = []
    sig_peak_array = []
    distance_modulus_array = []
    n_obs_peak_array = []
    n_obs_half_peak_array = []
    n_model_peak_array = []

    x_peak_array, y_peak_array, angsep_peak_array = region.find_peaks(data, distance_modulus)

    #for x_peak, y_peak, angsep_peak in itertools.izip(x_peak_array, y_peak_array, angsep_peak_array): [Python 2 implementation]
    for x_peak, y_peak, angsep_peak in zip(x_peak_array, y_peak_array, angsep_peak_array):
        # Aperture fitting
        print('Fitting aperture to hotspot...')
        ra_peaks, dec_peaks, r_peaks, sig_peaks, distance_moduli, n_obs_peaks, n_obs_half_peaks, n_model_peaks = region.fit_aperture(data, distance_modulus, x_peak, y_peak, angsep_peak)
        
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
    parser.add_argument('--config',type=str,required=False,default='defaults.yaml',
                        help='config file')
    parser.add_argument('--outfile',type=str,required=False,
                        help='Output file')
    parser.add_argument('--ra',type=float,required=True,
                        help='Target RA')
    parser.add_argument('--dec',type=float,required=True,
                        help='Target DEC')
    args = vars(parser.parse_args())

    with open(args['config'], 'r') as ymlfile:
        cfg = yaml.load(ymlfile)
        survey = simple.survey.Survey(cfg)

    #--------------------------------------------------------------------------

    region = simple.survey.Region(survey, args['ra'], args['dec'])
    print('Search coordinates: (RA, Dec) = ({:0.2f}, {:0.2f})'.format(region.ra, region.dec))
    print('Search healpixel: {} (nside = {})'.format(region.pix_center, region.nside))
    print('Healpixels: {}'.format(region.pix_neighbors))

    #--------------------------------------------------------------------------

    region.load_data()
    print('Found {} objects'.format(len(region.data)))
    if (len(region.data) == 0):
        print('Ending search.')
        exit()

    #--------------------------------------------------------------------------

    distance_modulus_search_array = np.arange(survey.moduli['start'], survey.moduli['end']+1e-10, survey.moduli['step'])

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
        ra_peaks, dec_peaks, r_peaks, sig_peaks, dist_moduli, n_obs_peaks, n_obs_half_peaks, n_model_peaks = search_by_distance(survey, region, distance_modulus)
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
        print('{:0.2f} sigma; (RA, Dec, d) = ({:0.2f}, {:0.2f}); r = {:0.3f} deg; d = {:0.1f}, mu = {:0.2f} mag), mc_source_id: {:0.2f}'.format(sig_peak_array[ii], 
                         ra_peak_array[ii], 
                         dec_peak_array[ii], 
                         r_peak_array[ii],
                         ugali.utils.projector.distanceModulusToDistance(distance_modulus_array[ii]),
                         distance_modulus_array[ii],
                         mc_source_id_array[ii]))
    
    # Write output
    if (len(sig_peak_array) > 0):
        write_output(survey.output['results_dir'], survey.catalog['nside'], region.pix_center, ra_peak_array, dec_peak_array,
                     r_peak_array, distance_modulus_array, 
                     n_obs_peak_array, n_obs_half_peak_array, n_model_peak_array, 
                     sig_peak_array, mc_source_id_array, 0, args['outfile'])
    else:
        print('No significant hotspots found.')
        nan_array = [np.nan]
        write_output(survey.output['results_dir'], survey.catalog['nside'], region.pix_center,
                     nan_array, nan_array, nan_array, nan_array, 
                     nan_array, nan_array, nan_array, nan_array,
                     [mc_source_id], 0, args['outfile'])
