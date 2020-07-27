#!/usr/bin/env python
"""
Perform simple binning search
"""
__author__ = "Sidney Mau"

import os
import sys
import glob
import argparse
import subprocess
import yaml
import copy
import time

import healpy as hp
import ugali.utils.healpix

import simple.search

#------------------------------------------------------------------------------

def submit_job(scratch_dir, cfg, infile):
    # Read filename style and center healpix from infile
    splt = infile.split('.fits')[0].split('_')
    pix = int(splt[-1])
    fbase = '_'.join(splt[:-1])+'_{:0>5n}.fits'

    # Initialize run directory
    rundir = "{}/rundir-{:0>5n}".format(scratch_dir, pix)
    subprocess.call("mkdir {}".format(rundir).split())

    # Copy over datafiles and grid selections
    subprocess.call("cp {} {}".format(infile, rundir).split())
    for p in hp.get_all_neighbours(cfg['catalog']['nside'], pix):
        fname = fbase.format(p)
        subprocess.call("cp {} {}".format(fname, rundir).split())
    grid_selection = glob.glob("{}/*{}.npz".format(cfg['grid']['grid_dir'], pix))[0]
    subprocess.call("cp {} {}".format(grid_selection, rundir).split())

    ra, dec = ugali.utils.healpix.pixToAng(cfg['catalog']['nside'], pix)
    ra, dec = round(ra, 2), round(dec, 2)
    subprocess.call("{}/submitscr.sh {:0>5n} {} {}".format(scratch_dir, pix, ra, dec).split())


#------------------------------------------------------------------------------

if __name__ == '__main__':
    # Construct argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument('--config',type=str,required=True,
                        help='config file')
    parser.add_argument('--njobs',type=int,required=False,
                        help='maximum jobs in queue')
    parser.add_argument('--gb',type=float,required=False,
                        help='maximum size (GB) of scratch run diretory')
    args = vars(parser.parse_args())

    if args['njobs'] is None and args['gb'] is None:
        cont = input("WARNING: No maximum job number or directory size specified. Are you sure you want to coninue? (y/n)")
        if ('y' not in cont) and ('Y' not in cont):
           sys.exit(0)


    with open(args['config'], 'r') as ymlfile:
        cfg = yaml.load(ymlfile)

    #--------------------------------------------------------------------------

    scratch_dir = "/scratch/mcnanna/simple_run"
    subprocess.call("mkdir {}".format(scratch_dir).split())
    
    # Write submit file in this directory
    with open('submitscr.sh') as f:
        s = f.read()
    lines = s.splitlines()
    for j in range(len(lines)):
        if 'transfer_input_files' in lines[j].lower():
            for module in cfg['catalog']['other'].split('&&'):
                lines[j] += ',{}.py'.format(module.strip())
    with open('{}/submitscr.sh'.format(scratch_dir), 'w') as f:
        f.write('\n'.join(lines))
    subprocess.call("chmod +x {}/submitscr.sh".format(scratch_dir).split())

    # Write init_py in this directory
    with open('init_python.sh') as f:
        p = f.read()
    with open('{}/init_python.sh'.format(scratch_dir), 'w') as f:
        f.write(p)
    subprocess.call("chmod +x {}/init_python.sh".format(scratch_dir).split())

    # Copy over search.py and fracdet file
    subprocess.call("cp {}/search.py {}".format(os.path.dirname(simple.search.__file__), scratch_dir).split())
    if cfg['survey']['fracdet'] is not None:
        subprocess.call("cp {} {}/fracdet.fits".format(cfg['survey']['fracdet'], scratch_dir).split())

    # Write new config file w/ modified directories. Copy any additional modules in the process
    new_cfg = copy.deepcopy(cfg)
    new_cfg['catalog']['dirname'] = '.'
    try:
        new_cfg['grid']['grid_dir'] = '.'
    except:
        pass
    if new_cfg['survey']['fracdet'] is not None:
        new_cfg['survey']['fracdet'] = '../fracdet.fits'
    new_cfg['output']['results_dir'] = '.'
    if new_cfg['catalog']['other'] is not None:
        modules = []
        for module in new_cfg['catalog']['other'].split('&&'):
            subprocess.call("cp {}.py {}".format(module.strip(), scratch_dir).split())
            modules.append('../{}'.format(module.strip()))
        new_cfg['catalog']['other'] = ' && '.join(modules)

    with open('{}/config.yaml'.format(scratch_dir), 'w') as newymlfile:
        yaml.dump(new_cfg, newymlfile)


    def clean_up():
        finished_dirs = ['/'.join(result.split('/')[:-1]) for result in glob.glob('{}/rundir*/outfile'.format(scratch_dir))]
        if len(finished_dirs) > 0:
            print('Cleaning up results for {} pix.'.format(len(finished_dirs)))
            for dirname in finished_dirs:
                subprocess.call("mv {} {}".format(dirname+'/results.npy', cfg['output']['results_dir']+'/results_nside_{}_{}.npy'.format(cfg['catalog']['nside'], dirname.split('-')[-1])).split())
                #subprocess.call("mv {} {}".format(dirname+'/outfile', cfg['output']['log_dir']+'/log_nside_{}_{}.log'.format(cfg['catalog']['nside'], dirname.split('-')[-1])).split())
                hostinfo = open('{}/hostinfo'.format(dirname), 'r')
                outfile = open('{}/outfile'.format(dirname), 'r')
                logfile = open(cfg['output']['log_dir']+'/log_nside_{}_{}.log'.format(cfg['catalog']['nside'], dirname.split('-')[-1]), 'w')
                logfile.write(hostinfo.read())
                logfile.write(outfile.read())
                hostinfo.close()
                outfile.close()
                logfile.close()
                subprocess.call("rm -r {}".format(dirname).split())
        return len(finished_dirs)

    # Read infiles, submit job for each
    infiles = glob.glob('{}/*.fits'.format(cfg['catalog']['dirname']))[::-1]
    # TODO: This list is reversed to run in tandem with a local session running forwards

    # Read already completed files based on results_dir, skip those
    donefiles = glob.glob('{}/*.npy'.format(cfg['output']['results_dir']))
    done_pix = [f.split('.npy')[0].split('_')[-1] for f in donefiles]

    submitted = 0
    finished = 0
    i = 0
    while i < len(infiles):
        pix = infiles[i].split('.fits')[0].split('_')[-1]
        if pix in done_pix:
            print("Skipping pix={}: already done".format(pix))
            i += 1
            continue
        
        # Check if space is used up or queue is full
        if args['gb'] is None:
            space_full = False
        else:
            out = str(subprocess.check_output("du -sh /scratch/mcnanna", shell=True))
            mem = out.split('\\t')[0].split("'")[-1]
            space_full = (mem[-1] == 'G' and float(mem[:-1]) >= args['gb'])
        if args['njobs'] is None:
            queue_full = False
        else:
            queue_full = (submitted - finished >= args['njobs'])

        if (not queue_full) and (not space_full):
            submit_job(scratch_dir, cfg, infiles[i])
            submitted += 1
            i += 1
            print('Sumbitted pix {} ({}/{})'.format(pix, submitted, len(infiles)-len(done_pix)))
            # Check for any completed jobs and copy the files to their permanent location:
            finished += clean_up()
        else: 
            # Wait until for more room in queue/memory
            if space_full:
                print('Space limit reached. Waiting for more room...')
            if queue_full:
                print('Job limit reached. Waiting for space in the queue...')
            while len(glob.glob('{}/rundir*/outfile'.format(scratch_dir))) == 0:
                time.sleep(60)
            # Copy results into permanent location, delete rundirs
            finished += clean_up()

    # Copy over remaining results as the jobs finish
    print('All jobs sumbitted. Waiting to clean them up.')
    while finished < submitted:
        time.sleep(60*10)
        finished += clean_up()
