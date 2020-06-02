#!/usr/bin/env python
"""
Perform simple binning search
"""
__author__ = "Sidney Mau"

import os
import glob
import argparse
import subprocess
import yaml

import ugali.utils.healpix

import simple.search

#------------------------------------------------------------------------------

def submit_job(cfgfile, cfg, ra, dec, pix):
    outfile = '{}/results_nside_{}_{}.npy'.format(cfg['output']['results_dir'], cfg['catalog']['nside'], pix)
    logfile = '{}/results_nside_{}_{}.log'.format(cfg['output']['log_dir'], cfg['catalog']['nside'], pix)
    batch = 'csub -n {} -o {} '.format(cfg['batch']['max_jobs'], logfile)
    command = 'python {}/search.py --config {} --ra {:0.2f} --dec {:0.2f} --outfile {}'.format(os.path.dirname(simple.search.__file__), cfgfile, ra, dec, outfile)
    command_queue = batch + command

    print(command_queue)
    subprocess.call(command_queue.split(' '), shell=False)

    return

#------------------------------------------------------------------------------

if __name__ == '__main__':
    # Construct argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument('--config',type=str,required=False,default='config.yaml',
                        help='config file')
    args = vars(parser.parse_args())

    with open(args['config'], 'r') as ymlfile:
        cfg = yaml.load(ymlfile)

    #--------------------------------------------------------------------------

    infiles = glob.glob ('{}/*.fits'.format(cfg['catalog']['dirname']))
    
    print('Pixelizing...')
    pix_nside = [] # Equatorial coordinates, RING ordering scheme
    for infile in infiles:
        pix_nside.append(int(infile.split('.fits')[0].split('_')[-1]))
    
    for ii in range(0, len(pix_nside)):
        ra, dec = ugali.utils.healpix.pixToAng(cfg['catalog']['nside'], pix_nside[ii])
    
        submit_job(args['config'], cfg, ra, dec, pix_nside[ii])
        print('({}/{})').format(ii+1, len(pix_nside))
