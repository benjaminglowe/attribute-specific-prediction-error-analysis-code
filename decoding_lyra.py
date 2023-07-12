# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 16:02:53 2022

@author: Benjamin G. Lowe (ben.lowe@mq.edu.au)

Script run on QUT's HPC cluster's (collectively named 'Lyra'). This is not
meant to be run directly by the user, but rather called with qsub_all.sh.
"""
#%% Importing important libraries
import os
import sys
import glob
import numpy as np
import pandas as pd
from scipy import io
from project_funcs import process_participant
from concurrent.futures import ProcessPoolExecutor, as_completed

#%% Defining important variables
root = ''                                                                       # BIDS root directory
data_path = f'{root}//derivatives//epo//MVPA'
output_path = f'{root}//derivatives//results//MVPA'
conditions = ['col_vio', 'siz_vio', 'ori_vio', 'control']
perms = 100                                                                     # number of permutations per participant
C = 1                                                                           # training cost parameter
k = 5                                                                           # k-folds
min_trials = 64                                                                 # number of trials needed to make it through rejection for participant to be analysed     
time_points = 181                                                               # number of time points in epoch
baseline = (-0.1, 0)                                                            # baseline window
os.chdir(root)                                                                  # changing directory to BIDS root
ch_picks = list(pd.read_csv('code//EEG_channels_64.csv')['labels'])             # channels to include in analysis
os.chdir(data_path)

#%% User inputs
fif_file = sys.argv[1]                                                          # fif file to analyse (given by bash script)

#%% Decoding participants
sub = 'sub-{}'.format(fif_file.split('-')[1])                                   # subject ID 
os.chdir(f'{root}//{sub}//beh')
beh_data = pd.read_csv(glob.glob('*.csv')[0])                                   # loading behavioural data into a dataframe
beh_data.dropna(axis=0, inplace=True, subset=['im_0_file'])                     # removing breaks from behavioural dataframe (nan rows)
beh_data.reset_index(inplace=True)                                              # reseting index for below loop
catch_indicies = []                                                             # list of catch trial indicies within behavioural dataframe
for i in range(len(beh_data)):
    if '(CATCH)' in beh_data.loc[i]['trial_label']:
        catch_indicies.append(i)
beh_data.drop(axis=0, index=catch_indicies, inplace=True)                       # removing catch trials from behavioural dataframe
beh_data.drop('index', axis=1, inplace=True)                                    # removing addition 'index' column from behavioural dataframe
os.chdir(data_path)
if __name__ == '__main__':
    vio_vs_cont = np.zeros((perms, 3, time_points, time_points))                # array for storing performance for each attribute violation vs. control trial subsets
    vio_vs_cont_SL = np.zeros((perms, 2, time_points, len(ch_picks)))           # array for storing weights corresponding with above models
    siz_vs_ori = np.zeros((perms, time_points, time_points))
    siz_vs_ori_SL = np.zeros((perms, time_points, len(ch_picks)))
    TL = np.zeros((perms, 2, 2, time_points, time_points))                      # array for storing transfer learning matricies
    main_effect = np.zeros((perms, time_points, time_points))
    main_effect_SL = np.zeros((perms, time_points, len(ch_picks)))
    with ProcessPoolExecutor() as executor:                                     # initialising parallel fork
        results = [executor.submit(process_participant, fif_file, sub, 
                                   conditions, baseline, ch_picks, 
                                   min_trials, beh_data, C, 
                                   k, rep_avg=False) for i in range(perms)]     # running decoding pipeline in parallel across cores
        for i, future in enumerate(as_completed(results)):
            try:                                                                # try storing the results
                result = future.result()
                vio_vs_cont[i, :, :, :] = result[0]
                vio_vs_cont_SL[i, :, :, :] = result[1]
                siz_vs_ori[i, :, :] = result[2]
                siz_vs_ori_SL[i, :, :] = result[3]
                TL[i, :, :, :, :] = result[4]
                main_effect[i, :, :] = result[5]
                main_effect_SL[i, :, :] = result[6]
                save_data = True
            except:                                                             # if that fails, part the participant as bad
                save_data = False
            print(f'save_data = {save_data}') 
    if save_data:
        os.chdir(output_path)
        mat_file = {'vio_vs_cont': vio_vs_cont,
                    'vio_vs_cont_SL': vio_vs_cont_SL,
                    'siz_vs_ori': siz_vs_ori,
                    'siz_vs_ori_SL': siz_vs_ori_SL,
                    'TL': TL,
                    'main_effect': main_effect,
                    'main_effect_SL': main_effect_SL}
        io.savemat(f'{sub}-decode.mat', mat_file)
        