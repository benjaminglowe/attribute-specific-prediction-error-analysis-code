# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 12:13:35 2021

@author: Benjamin G. Lowe (ben.lowe@mq.edu.au)

Preprocessing script. Assumes that data are formatted according to BIDS. Set
MVPA_data to True for decoding analysis.

This can be run locally, but it was run on QUT's HPC clusters (collectively
called 'Lyra') by running qsub_preprocess_all.sh.
"""
#%% Importing important libraries
import os
import mne
import glob
import numpy as np
import pandas as pd
from project_funcs import stayed_vigilant

#%% Defining important variables
path = ''                                                                       # path to BIDS directory
montage = mne.channels.make_standard_montage(kind='biosemi64')                  # EEG montage used during recording (important for interpolations)
ref_channels = 'average'                                                        # Re-referencing "electrode"
tshift = 1/60*2*-1                                                              # accounting for image/trigger asynchronoy
cut_off = 0.9                                                                   # percentage of catch trials a participant must have correctly responded
bad_participants = {}                                                           # empty dictionary for storing bad participant
os.chdir(path)
subs = glob.glob('sub*')
MVPA_data = False

#%% Dataset specific parameters
if MVPA_data:
    Fs = 200                                                                    # resampling rate
    filt_method = 'fir'                                                         # filtering method
    l_cut = 0.1                                                                 # high-pass filter
    h_cut = 100                                                                 # low-pass filter
    notch = None                                                                # notch filter frequencey
    tmin = -0.2                                                                 # epoch time prior to event
    tmax = 0.7                                                                  # epoch time post event
else:
    Fs = 200                                                                    # resampling rate
    filt_method = 'fir'                                                         # filtering method
    l_cut = 0.1                                                                 # high-pass filter
    h_cut = 100                                                                 # low-pass filter
    notch = 50                                                                  # notch filter frequencey
    tmin = -0.2                                                                 # epoch time prior to event
    tmax = 0.7                                                                  # epoch time post event

#%% Defining bad channels per participant
bad_channels = {'sub-01': ['F3'],
                'sub-02': [],
                'sub-03': [],
                'sub-04': [],
                'sub-05': [],
                'sub-06': [],
                'sub-07': ['T7'],
                'sub-08': ['TP7', 'P9'],
                'sub-09': [],
                'sub-10': [],
                'sub-11': [],
                'sub-12': ['P5'],
                'sub-13': ['TP7'],
                'sub-14': [],
                'sub-15': [],
                'sub-16': [],
                'sub-17': [],
                'sub-18': [],
                'sub-19': [],
                'sub-20': [],
                'sub-21': [],
                'sub-22': ['F7', 'FT7', 'T7'],
                'sub-23': ['FT7', 'TP8', 'PO8'],
                'sub-24': ['C6'],
                'sub-25': ['PO4'],
                'sub-26': ['T7', 'Fpz'],
                'sub-27': [],
                'sub-28': ['F8', 'AF8', 'AF7'],
                'sub-29': [],
                'sub-30': ['P10', 'PO8', 'P8'],
                'sub-31': [],
                'sub-32': ['TP8', 'T8', 'Fp2'],
                'sub-33': ['F4'],
                'sub-34': [],
                'sub-35': [],
                'sub-36': [],
                'sub-37': ['P6', 'Fp1', 'AF7']}

#%% Loading in participant data
for sub in subs:
    if sub == 'sub-11':
        continue
    raw = mne.io.read_raw_brainvision(f'{sub}//eeg//{sub}_task-CTP_eeg.vhdr',   # loading in participant data
                                      preload=True)
    beh = pd.read_csv(f'{sub}//beh//{sub}_task-CTP_beh.csv')                    # loading in behavioural data
    if stayed_vigilant(beh, cut_off) == False:                                  # if participant failed to pay attention to task
        bad_participants[sub] = 'not vigilant'
        continue

#%% Interpolating
    raw.set_montage(montage=montage)
    raw.info['bads'] = bad_channels[sub]
    if len(raw.info['bads']) > 4:
        bad_participants[sub] = 'too many bad channels'
        continue
    raw.interpolate_bads(reset_bads=True)

#%% Re-reference
    raw.set_eeg_reference(ref_channels=ref_channels)

#%% Filtering
    if l_cut != None or h_cut != None:
        raw.filter(l_cut, h_cut, method=filt_method)
    if notch != None:
        raw.notch_filter(freqs=notch, method=filt_method)

#%% Resampling
    if Fs != None:
        raw.resample(Fs)

#%% Extracting events
    conditions = {'col_vio': 'Stimulus/S   10',                                 # condition to event map
                  'siz_vio': 'Stimulus/S   20',
                  'ori_vio': 'Stimulus/S   30',
                  'control': 'Stimulus/S   40'}
    events, event_id = mne.events_from_annotations(raw)
    for e_key in event_id.keys():
        for c_key in conditions.keys():
            if e_key == conditions[c_key]:
                conditions[c_key] = event_id[e_key]

#%% Epoching
    epochs = mne.Epochs(raw, events=events, event_id=conditions, preload=True,
                        reject=None, event_repeated='drop', baseline=None,
                        tmin=tmin-tshift,
                        tmax=tmax-tshift)
    epochs.shift_time(tshift=tshift, relative=True)

#%% Saving data
    if MVPA_data:
        epochs.save(f'derivatives//epo//MVPA//{sub}-epo.fif', overwrite=True)
    else:
        epochs.save(f'derivatives//epo//univariate//{sub}-epo.fif',
                    overwrite=True)
