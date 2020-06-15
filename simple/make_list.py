#!/usr/bin/env python
"""
Compile candidate list from results_dir
"""
__author__ = "Sidney Mau"

import glob
import yaml
import argparse

from astropy.io import fits
import numpy as np
import csv
import fitsio

#------------------------------------------------------------------------------

if __name__ == '__main__':
    # Construct argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument('--config',type=str,required=True,help='config file')
    parser.add_argument('--outname',type=str,required=False,help='Output .fits file name',default='candidate_list.fits')
    args = vars(parser.parse_args())

    with open(args['config'], 'r') as ymlfile:
        cfg = yaml.load(ymlfile)

    #--------------------------------------------------------------------------

    # Concatenate results from results_dir into one array
    results = []
    for file in glob.glob('{}/*.npy'.format(cfg['output']['results_dir'])):
        result = np.load(file)
        for row in result:
            results.append(row)
    data = np.asarray(results)

    # Write fits output
    outname = args['outname']
    if '.' not in outname:
        outname += '.fits'
    fits.writeto(outname, data, overwrite=True)
    
    # Diagnostic output
    data = fitsio.read(outname)
    print("{} hotspots found.".format(len(data)))
    cut_0 = (data['SIG'] > 5.5)
    print("{} hotspots found with SIG > 5.5.".format(len(data[cut_0])))
    cut_1 = (data['SIG'] > 10)
    print("{} hotspots found with SIG > 10.".format(len(data[cut_1])))
    cut_2 = (data['SIG'] > 15)
    print("{} hotspots found with SIG > 15.".format(len(data[cut_2])))
    cut_3 = (data['SIG'] > 20)
    print("{} hotspots found with SIG > 20.".format(len(data[cut_3])))
    cut_4 = (data['SIG'] > 25)
    print("{} hotspots found with SIG > 25.".format(len(data[cut_4])))
    cut_5 = (data['SIG'] >= 37.5)
    print("{} hotspots found with SIG >= 37.55".format(len(data[cut_5])))
