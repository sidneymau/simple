#!/usr/bin/env python
"""
"""
__author__ = "Sidney Mau"

# Python libraries
import os
import yaml
import argparse

# Simple Libraries
import simple.survey

#----------------------------------------------------------

def init_dirs(survey):
    results_dir = os.path.join(os.getcwd(), survey.output['results_dir'])
    if not os.path.exists(results_dir):
        os.mkdir(results_dir)

    log_dir = os.path.join(os.getcwd(), survey.output['log_dir'])
    if not os.path.exists(log_dir):
        os.mkdir(log_dir)

    save_dir = os.path.join(os.getcwd(), survey.output['save_dir'])
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)

if __name__ == '__main__':
    ## Construct argument parser
    #parser = argparse.ArgumentParser()
    #parser.add_argument('-i','--infile',type=str,required=True,
    #                    help='Input file')
    #args = vars(parser.parse_args())
    #args['infile']

    with open('defaults.yaml', 'r') as ymlfile:
        cfg = yaml.load(ymlfile)
        survey = simple.survey.Survey(cfg)

    init_dirs(survey)

    import pdb;pdb.set_trace()
