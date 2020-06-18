#!/usr/bin/env python

import numpy as np
import yaml
import subprocess
import argparse
import glob

import ugali.utils.healpix
import ugali.utils.projector

parser = argparse.ArgumentParser()
parser.add_argument('--config', type=str, required=True)
args = vars(parser.parse_args())

with open(args['config'], 'r') as ymlfile:
    cfg = yaml.load(ymlfile)

subprocess.call('mkdir -p {}'.format(cfg['grid']['grid_dir']).split())
infiles = glob.glob('{}/*.fits'.format(cfg['catalog']['dirname']))

pix_nside = [] # Equatorial coordinates, RING ordering scheme
for infile in infiles:
    pix_nside.append(int(infile.split('.fits')[0].split('_')[-1]))

area = cfg['grid']['delta_x']**2
bins = np.arange(-cfg['grid']['bin_edge'], cfg['grid']['bin_edge']+1.e-10, cfg['grid']['delta_x'])
centers = 0.5 * (bins[0: -1] + bins[1:])

yy, xx = np.meshgrid(centers, centers)
xxflat, yyflat = xx.flatten(), yy.flatten()

for pix_nside_select in pix_nside:
    ra, dec = ugali.utils.healpix.pixToAng(cfg['catalog']['nside'], pix_nside_select)
    proj = ugali.utils.projector.Projector(ra, dec)

    rara, decdec = proj.imageToSphere(xxflat, yyflat)
    cutcut = (ugali.utils.healpix.angToPix(cfg['catalog']['nside'], rara, decdec) == pix_nside_select).reshape(xx.shape)
    np.savez_compressed('{}/grid_selection_{}'.format(cfg['grid']['grid_dir'], pix_nside_select), cutcut)
    

