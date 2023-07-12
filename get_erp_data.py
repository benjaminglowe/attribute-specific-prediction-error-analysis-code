#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 16:30:02 2023

@author: Benjamin G. Lowe (ben.lowe@mq.edu.au)

Script to be run by user for getting ERP data. Assumes that preprocessing.py
has been run with MVPA_data set to false.

This script was run locally.
"""
#%% Importing important libraries
import os
import glob
import mne
import numpy as np
from scipy import io

#%% Defining important variables
root = ''                                                                       # root of BIDS directory
epo_dir = f'{root}/derivatives/epo/univariate'
fif_files = glob.glob(f'{epo_dir}/*.fif')
fif_files.sort()
N = len(fif_files)
results_dir = f'{root}/derivatives/results/ERPs'
conditions = ['col_vio', 'siz_vio', 'ori_vio', 'control']
baseline = (-0.1, 0)
time_points = 181
n_chans = 64

#%% Making requiring directories
if not os.path.isdir(results_dir):
    os.mkdir(results_dir)

#%%
evokeds = np.zeros((N, len(conditions), n_chans, time_points))                  # empty array for storing participant-level, evoked data
for i, fif_file in enumerate(fif_files):
    epochs = mne.read_epochs(fif_file)                                          # reading in the data
    epochs.apply_baseline(baseline=baseline)                                    # applying a baseline correction
    epochs.drop_channels('Status')

    for j, condition in enumerate(conditions):
        evokeds[i, j, :, :] = epochs[condition].average().get_data()

#%% Organising data
mat_data = {}
for i, condition in enumerate(conditions):
    mat_data[condition] = evokeds[:, i, :, :]

#%% Saving data
io.savemat(f'{results_dir}/erps.mat', mat_data)
