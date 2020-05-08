#!/usr/bin/env python
"""
"""
__author__ = "Sidney Mau"

# Python libraries
import os
import yaml

# Simple Libraries
import simple.survey

#----------------------------------------------------------

def init_dir(survey):
    results_dir = os.path.join(os.getcwd(), survey['output']['results_dir'])
    if not os.path.exists(results_dir):
        os.mkdir(results_dir)

    log_dir = os.path.join(os.getcwd(), survey['output']['log_dir'])
    if not os.path.exists(log_dir):
        os.mkdir(log_dir)

    save_dir = os.path.join(os.getcwd(), survey['output']['save_dir'])
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)

if __name__ == '__main__':
    with open('defaults.yaml', 'r') as ymlfile:
        cfg = yaml.load(ymlfile)
        survey = simple.survey.Survey(cfg)

    import pdb;pdb.set_trace()
