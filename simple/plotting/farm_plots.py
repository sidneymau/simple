#!/usr/bin/env python
"""
Create simple binner style plots for ugali or simple candidate lists
"""
__author__ = "Sidney Mau"

import os
import subprocess
import numpy as np
import argparse
import yaml

import simple.survey
import simple.plotting.diagnostic_plots

############################################################

parser = argparse.ArgumentParser()
parser.add_argument('--config', type=str, required=True, help='config file')
parser.add_argument('--infile', type=str, required=True, help='candidate list to plot')
parser.add_argument('--outdir', type=str, required=False, help='directory for plots', default='plot_dir')
parser.add_argument('--sig_cut',type=str, required=False, help='significance cut', default=5.5)
parser.add_argument('--jobs'   ,type=int, required=False, help='number of simultaneous jobs', default=20)
args = vars(parser.parse_args())

with open(args['config'], 'r') as ymlfile:
    cfg = yaml.load(ymlfile)
    survey = simple.survey.Survey(cfg)

save_dir = os.path.join(os.getcwd(), args['outdir'])
if not os.path.exists(save_dir):
    os.mkdir(save_dir)

log_dir = os.path.join(save_dir, 'log_dir')
if not os.path.exists(log_dir):
    os.mkdir(log_dir)


print('Plotting hotspots with sig > {}'.format(args['sig_cut']))

candidate_list = np.load(args['infile'])
try: # simple
    candidate_list = candidate_list[candidate_list['SIG'] > args['sig_cut']]
except: # ugali
    candidate_list = candidate_list[candidate_list['TS'] > 25]

print('{} candidates found...').format(len(candidate_list))

############################################################

#for candidate in [candidate_list[:10]]:
for candidate in candidate_list:
    keys = ['sig', 'ra', 'dec', 'mod', 'r', 'n_obs', 'n_model']
    params = dict.fromkeys(keys)
    if candidate is not None:
        try: # simple
            params['sig'] = round(candidate['SIG'], 2)
        except: # ugali
            params['sig'] = round(candidate['TS'], 2)
        params['ra']      = round(candidate[survey.catalog['basis_1']], 2)
        params['dec']     = round(candidate[survey.catalog['basis_2']], 2)
        params['mod']     = round(candidate['MODULUS'], 2)
        params['r']       = round(candidate['R'], 3)
        params['n_obs']   = candidate['N_OBS']
        params['n_model'] = candidate['N_MODEL']
    sig, ra, dec, mod, r, n_obs, n_model = [params[key] for key in keys]
    
    logfile = '{}/candidate_{}_{}.log'.format(log_dir, ra, dec)
    #batch = 'csub -n {} -o {} '.format(jobs, logfile)
    batch = 'csub -n {} -o {} --host all '.format(args['jobs'], logfile) # testing condor updates
    abspath = os.path.abspath(simple.plotting.diagnostic_plots.__file__)
    command = 'python {} --config {} --outdir {} --sig {:0.2f} --ra {:0.2f} --dec {:0.2f} --r {:0.3f} --modulus {:0.2f} --n_obs {:0.2f} --n_model {:0.2f}'.format(
            abspath, args['config'], args['outdir'], sig, ra, dec, r, mod, n_obs, n_model)
    command_queue = batch + command

    print(command_queue)
    subprocess.call(command_queue.split(' '), shell=False)
