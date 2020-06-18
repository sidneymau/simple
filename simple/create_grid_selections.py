#!/usr/bin/env python

import numpy as np
import yaml
import subprocess
import argparse
import glob

import ugali.utils.healpix
import ugali.utils.projector

import simple.survey

parser = argparse.ArgumentParser()
parser.add_argument('--config', type=str, required=True)
args = vars(parser.parse_args())

with open(args['config'], 'r') as ymlfile:
    cfg = yaml.load(ymlfile)
    survey = simple.survey.Survey(cfg)

subprocess.call('mkdir -p {}'.format(survey.grid['grid_dir']).split())
infiles = glob.glob('{}/*.fits'.format(survey.catalog['dirname']))

pix_nside = [] # Equatorial coordinates, RING ordering scheme
for infile in infiles:
    pix_nside.append(int(infile.split('.fits')[0].split('_')[-1]))

area = survey.grid['delta_x']**2
bins = np.arange(-survey.grid['bin_edge'], survey.grid['bin_edge']+1.e-10, survey.grid['delta_x'])
centers = 0.5 * (bins[0: -1] + bins[1:])

yy, xx = np.meshgrid(centers, centers)
xxflat, yyflat = xx.flatten(), yy.flatten()

for pix_nside_select in pix_nside:
    ra, dec = ugali.utils.healpix.pixToAng(survey.catalog['nside'], pix_nside_select)
    proj = ugali.utils.projector.Projector(ra, dec)

    rara, decdec = proj.imageToSphere(xxflat, yyflat)
    cutcut = (ugali.utils.healpix.angToPix(survey.catalog['nside'], rara, decdec) == pix_nside_select).reshape(xx.shape)
    np.savez_compressed('{}/grid_selection_{}'.format(survey.grid['grid_dir'], pix_nside_select), cutcut)


