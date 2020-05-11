#!/usr/bin/env python
"""
Perform simple binning search
"""
__author__ = "Sidney Mau"

import glob
import subprocess
import yaml

import ugali.utils.healpix

#------------------------------------------------------------------------------

def submit_job(cfg, ra, dec, pix):
    outfile = '{}/results_nside_{}_{}.txt'.format(results_dir, nside, pix)
    logfile = '{}/results_nside_{}_{}.log'.format(log_dir, nside, pix)
    batch = 'csub -n {} -o {} '.format(jobs, logfile)
    command = 'python {}/search.py --config {} --ra {:0.2f} --dec {:0.2f} --outfile {} --logfile {}'.format(simple_dir, cfgfile, ra, dec, mc_source_id, outfile, logfile)
    command_queue = batch + command

    print(command_queue)
    subprocess.call(command_queue.split(' '), shell=False)

    return

#------------------------------------------------------------------------------

if __name__ == '__main__':
    # Construct argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument('--config',type=str,required=False,default='defaults.yaml',
                        help='config file')
    args = vars(parser.parse_args())

    with open(args['config'], 'r') as ymlfile:
        cfg = yaml.load(ymlfile)

    #--------------------------------------------------------------------------

    infiles = glob.glob ('{}/*.fits'.format(cfg['datadir']))
    
    print('Pixelizing...')
    pix_nside = [] # Equatorial coordinates, RING ordering scheme
    for infile in infiles:
        pix_nside.append(int(infile.split('.fits')[0].split('_')[-1]))
    
    for ii in range(0, len(pix_nside)):
        ra, dec = ugali.utils.healpix.pixToAng(cfg['nside'], pix_nside[ii])
    
        submit_job(args['config'], ra, dec, pix_nside[ii])
        print('({}/{})').format(ii+1, len(pix_nside))
