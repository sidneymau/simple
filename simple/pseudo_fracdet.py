#!/usr/bin/env python
"""
Pseudo-fracdet map
"""
__author__ = "Sidney Mau"

import argparse
import os
import glob
import yaml
import numpy as np
import healpy as hp
import fitsio as fits

parser = argparse.ArgumentParser()
parser.add_argument('--config', type=str, required=True)
args = vars(parser.parse_args())

with open(args['config'], 'r') as ymlfile:
    cfg = yaml.load(ymlfile)


############################################################

infiles = glob.glob ('{}/*.fits'.format(cfg['catalog']['dirname']))

nside = 2048
pix = []
for infile in infiles:
    print('loading {}'.format(infile))
    data = fits.read(infile, columns=[cfg['catalog']['basis_1'], cfg['catalog']['basis_2']])
    p = hp.ang2pix(nside, data[cfg['catalog']['basis_1']], data[cfg['catalog']['basis_2']], lonlat=True)
    pix.append(np.unique(p))

print('Constructing map')
pix = np.concatenate(pix)
pix = np.unique(pix)
coverage_map = np.tile(hp.UNSEEN, hp.nside2npix(nside))
coverage_map[pix] = 1

print('Writing output')
result = '{}_pseudo_fracdet.fits.gz'.format(cfg['survey']['name'])
hp.write_map(result, coverage_map)
