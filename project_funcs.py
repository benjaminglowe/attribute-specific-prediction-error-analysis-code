# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 12:24:47 2021

@author: Benjamin G. Lowe (ben.lowe@mq.edu.au)

Functions written as part of Lowe's et al. decoding pipeline.
"""
#%% Importing important functions
import mne
import numpy as np
import pandas as pd
from math import floor
from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler as scaling_method

#%% Stayed vigilant
def stayed_vigilant(behav_data, cut_off):
    """
    Testing whether the participant was responded to a sufficient number of
    catch trials whilst viewing the paradigm.

    Parameters
    ----------
    behav_data : pandas DataFrame
        The participant's behavioural data output as read by pandas.
    cut_off : float
        Proportion of catch trials participant had to respond to in order to be
        considered vigilant.

    Returns
    -------
    vigilant : bool
        Statement as to whether or not the participant was deemed as vigilant.

    """
    catch_resps = []
    vigilant = False
    for i in range(len(behav_data)):
        if type(behav_data.loc[i]['trial_label']) != str:
            continue
        elif '(CATCH)' in behav_data.loc[i]['trial_label']:
            catch_resps.append(behav_data.loc[i]['key_resp_trial.corr'])
    if sum(catch_resps) >= len(catch_resps)*cut_off:
        vigilant = True

    return vigilant

#%% Renaming trial events to be more specific
def find_new_codes(epochs_events, behav_data):
    """
    Function for making trigger codes more specific based on:
        1. The trial condition
        2. A numerical representation of the final stimulus within a sequence

    This function was written so that we could later match for final stimulus
    representations between conditions. Thus, if a signal was decoded, it
    low-level differences in feature presentation were controlled for.

    Parameters
    ----------
    epochs_events : array
        A two-dimensional array specifying the timing of events. This should be
        called using the events attribute of an Epochs mne object.
    behav_data : pandas DataFrame
        The participant's behavioural data output as read by pandas.

    Returns
    -------
    epochs_events : array
        A two-dimensional array specifying the timing of events with more
        specific event triggers.

    """
    trial_labels = list(behav_data['trial_label'])                              # reading in specific trial labels from behavioural data sheet
    new_codes = []                                                              # list for storing new codes
    for i in range(len(trial_labels)):                                          #
        trial_label = trial_labels[i]                                           # defining the current label within iteration
        if trial_label.split(': ')[0] == 'col vio':                             # converting condition name to a number
            condition = '10'
        elif trial_label.split(': ')[0] == 'siz vio':
            condition = '20'
        elif trial_label.split(': ')[0] == 'ori vio':
            condition = '30'
        elif trial_label.split(': ')[0] == 'control':
            condition = '40'
        final_rep = trial_label.split(': ')[1].split('; ')                      # list describing the trajectory
        new_code = []                                                           # bin to store new code components in
        for entry in final_rep:
            new_code.append(''.join(entry.split('_')[1:3]))                     # removes anything that's not a number from string
        new_code = ''.join(new_code)                                            # joins components together
        new_codes.append(int(''.join([condition, new_code])))                   # joins together condition label and components for new code and appends them to the main list
    assert len(new_codes) == len(trial_labels) == len(epochs_events)            # makes sure no funny business is going on
    epochs_events[:, 2] = np.array(new_codes)                                   # changes column of event array to match newly created main trial list

    return epochs_events

#%% Finding shared trials between conditions
def find_unmatched_trials(epochs_events):
    """
    Function written to search for trials missing their physical representation
    match across all conditions. Output are epoch indicies missing a match.

    Parameters
    ----------
    epochs_events : MNE events array

    Returns
    -------
    drop_inds : list
        List of epoch indicies to drop.

    """
    col_vios = []                                                               # bins for storing specific condition events
    siz_vios = []
    ori_vios = []
    controls = []
    overlap = []                                                                # bin for overlap
    for i in range(len(epochs_events)):
        event = str(epochs_events[i, 2])                                        # turns current event into a string
        if event[:2] == '10':                                                   # appends appropriate bin with final rep if event belongs to that condition
            col_vios.append(event[2:len(event)])
        elif event[:2] == '20':
            siz_vios.append(event[2:len(event)])
        elif event[:2] == '30':
            ori_vios.append(event[2:len(event)])
        elif event[:2] == '40':
            controls.append(event[2:len(event)])

    col_vios = np.array(col_vios)                                               # converts lists into numpy arrays
    siz_vios = np.array(siz_vios)
    ori_vios = np.array(ori_vios)
    controls = np.array(controls)

    for trial in set(controls):                                                 # for every final rep (unique trial) that's in the control condition
        if trial in col_vios and trial in siz_vios and trial in ori_vios:       # if the final rep exists in all the control trials
            trial_counts = []                                                   # number of reps that trial has per condition
            trial_counts.append(len(np.where(col_vios == trial)[0]))
            trial_counts.append(len(np.where(siz_vios == trial)[0]))
            trial_counts.append(len(np.where(ori_vios == trial)[0]))
            trial_counts.append(len(np.where(controls == trial)[0]))
            for i in range(min(trial_counts)):                                  # minimum number of times this final rep has appeared per trial
                overlap.append(str(trial))                                      # depending on the number found within the line above, append this final rep to the overlap list that many times
    overlap.sort()                                                              # sorting the list into ascending order so that it is easier to work with

    epochs_to_keep = []
    for final_rep in overlap:                                                   # for every final rep within the overlap list, append epochs_to_keep with what that rep would look like as a trial code for each condition (after converting string to an integer)
        epochs_to_keep.append(int(''.join(['10', final_rep])))
        epochs_to_keep.append(int(''.join(['20', final_rep])))
        epochs_to_keep.append(int(''.join(['30', final_rep])))
        epochs_to_keep.append(int(''.join(['40', final_rep])))

    epochs_to_keep_ind = []                                                     # list of indicies of epochs to keep
    for epoch in epochs_to_keep:
        tmp = list(np.where(epochs_events[:, 2] == epoch)[0])                   # list of all trial indicies with that final rep
        ind = tmp[np.random.RandomState().randint(0, len(tmp))]                 # pulling a random index from the above list
        if ind not in epochs_to_keep_ind:                                       # store epoch indicies into list of epochs to keep (ensuring there's not double-ups)
            epochs_to_keep_ind.append(ind)
        else:
            while ind in epochs_to_keep_ind:                                    # resampling tmp until a unique index is called
                ind = tmp[np.random.RandomState().randint(0, len(tmp))]
            epochs_to_keep_ind.append(ind)
    epochs_to_keep_ind = set(epochs_to_keep_ind)                                #
    drop_inds = list(set(range(len(epochs_events))) - epochs_to_keep_ind)       # list of epoch indicies not in epochs_to_keep_ind
    assert len(epochs_to_keep_ind) + len(drop_inds) == len(epochs_events)       #

    return drop_inds

#%% Extracting main effect data
def get_main_effect_data(epochs, all_vio, control):
    """
    Function written to extract arrays of data used for the main effect
    comparison. Extracted data are stored into a dictionary. The matched
    indicies within the first dimension of each array correspond to data evoked
    by physically identical stimuli

    Parameters
    ----------
    epochs : MNE-Python epochs object

    Returns
    -------
    main_effect_data : dictionary

    """
    print('Extracting main effect data ...')

    epochs_events = epochs.events

    vio_events = list()                                                         # list for storing violaiton event types
    con_events = list()                                                         # list for storing control event types
    for i in range(epochs_events.shape[0]):
        tmp = list(str(epochs_events[i, 2]))                                    # list to index conditions and trial type -> first number indexes condition
        if tmp[0] != '4':                                                       # 4 indexes control trial
            tmp[0] = '1'                                                        # change to a 1 if not a 4. collapsing across violation types
            vio_events.append(''.join(tmp))
        else:
            con_events.append(''.join(tmp))

    assert len(vio_events) == len(con_events)                                   # making sure assumptions were met -> equal number of control and violation events
    vio_events = np.array(vio_events)
    con_events = np.array(con_events)

    vio_inds = list()                                                           # lists for storing trial indicies
    con_inds = list()
    labels = range(111112)                                                      # a range of trial labels to iterate through
    for label in labels:                                                        # iterating through range of labels

        exp_label = f'10{str(label).zfill(6)}'                                  # violation label
        vios = np.where(vio_events==exp_label)[0]                               # find when this occurred
        np.random.RandomState().shuffle(vios)                                   # randomise with every permutation
        vio_inds.extend(vios)                                                   # storing in above list

        exp_label = f'40{str(label).zfill(6)}'                                  # control label
        cons = np.where(con_events==exp_label)[0]                               # find when this occurred
        np.random.RandomState().shuffle(cons)                                   # randomise with every permutation
        con_inds.extend(cons)                                                   # storing in above list

    assert len(vio_inds) == len(vio_events)                                     # assumption checking
    assert len(vio_inds) == len(con_inds)
    assert all_vio.shape == control.shape

    main_effect_data = {'all_vio': all_vio[vio_inds, :, :],                     # storing ordered data into dictionary
                        'control': control[con_inds, :, :]}

    return main_effect_data

#%% Average across trial reps
def average_reps(epochs):
    """
    Function written to average across repetitions of the same trial. This is
    used to boost signal-to-noise during decoding.

    This transformation also re-orders epochs so that indicies are matched for
    physical characteristics across conditions. This is important for
    controlled cross-validation.

    Parameters
    ----------
    epochs : Epochs object within the MNE library

    Returns
    -------
    avg_data : dictionary
        Dictionary containing each condition's data after the transformation.
        Conditions (dictionary keys) are: col_vio, siz_vio, ori_vio, and
        control.
    """
    trial_types = list(set(epochs.events[:, 2]))                                # calls for each unique event type within the event list that made it through rejection
    trial_types.sort()
    chans, time_points = epochs.get_data().shape[1:]
    cond_data = np.zeros((len(trial_types), chans, time_points))
    for i, trial_type in enumerate(trial_types):
        type_inds = np.where(epochs.events == trial_type)[0]                    # finds indicies of particular trial type within epoch events
        cond_data[i, :, :] = epochs.get_data()[type_inds, :, :].mean(axis=0)    # averages across repeats of same trial type and stores within numpy array
    trials_n = len(trial_types)
    avg_data = {'col_vio': cond_data[:int(trials_n/4), :, :],
                'siz_vio': cond_data[int(trials_n/4):int((trials_n/4)*2), :, :],
                'ori_vio': cond_data[int((trials_n/4)*2):int((trials_n/4)*3), :, :],
                'control': cond_data[int((trials_n/4)*3):int((trials_n/4)*4), :, :]}

    return avg_data

#%% Reordering epochs within each condition
def find_reorder_indicies(epochs):
    """
    Finds a new order for epochs so that a given index within each condition
    has the same physical characteristics. This is important for controlled
    cross-validation.

    Parameters
    ----------
    epochs : Epochs object within MNE library

    Returns
    -------
    reordered_inds : list
        List of epoch indicies corresponding to the new order that allows for
        matching for physical stimulus presentation between conditions.
    """
    reordered_inds = []                                                         # "master list" for new epoch order
    labels = range(111112)                                                      # range of final (stimulus) rep labels (includes unused labels)
    condition_codes = [10, 20, 30, 40]
    for condition in condition_codes:                                           # iterating through conditions
        for label in labels:                                                    # iterating through range of labels
            label = str(label).zfill(6)                                         # padding out label
            condition_label = int('{}{}'.format(condition, label))              # joining condition label with final rep label
            inds = np.where(epochs.events == condition_label)[0]                # row indicies where condition label appears within events array
            for i in range(len(inds)):
                rand_i = np.random.RandomState().randint(0, len(inds))          # choosing a random index to call from inds
                reordered_inds.append(inds[rand_i])                             # appending above index to master list
                inds = np.delete(inds, rand_i)                                  # removing the randomly called index from inds

    return reordered_inds

#%% Decoding functions
def project_weights(X, clf):
    """
    Weight projection function (see Grootswager et al., 2017)

    Parameters
    ----------
    X : array
        Numpy array of data used to train the classifier. Rows correspond to
        observations, columns correspond to features.
    clf : sklearn classifier object
        Classifier fitted to training data. It is assumed that this is a linear
        classifier.

    Returns
    -------
    weights : array
        array corresponding to the projected weights (reconstructed activation
        pattern) from training data.
    """
    cov_X = np.cov(X, rowvar=False)
    w = clf.coef_.squeeze()
    A = np.dot(cov_X, w)

    return np.dot(A, 1/np.cov(np.dot(X, w)))

def linear_svm_decode(X, y, m, C=1, k=10):
    """
    Function written for linear decoding of evoked EEG data using a support-
    vector classifier.

    Parameters
    ----------
    X : array
        Data to decode (trials, channels, time points). It is assumed that
        there are an equal (balanced) number of trials per condition and that
        they are organised such that the first half of the data corresponds to
        one condition and vice versa.
    y : array or list
        Label corresponding to each trial.
    m : int
        Number of cases for each label.
    C : int, optional
        SVM cost parameter. The default is 1.
    k : int, optional
        Number of cross-validation folds. The default is 10.

    Returns
    -------
    time_series : array
        Decoding performance time series averaged across cross-validation 
        folds.
        
    weights : array
        Decoding weights time series averaged across cross-validation folds.

    """
    assert int(len(y)/len(list(set(y)))) == int(m), 'unequal number of cases per condition.'
    output = np.zeros((X.shape[2], k))                                          # time points by k-folds
    weights = np.zeros((X.shape[2], X.shape[1], k))
    fold_size = floor(m/k)                                                      # number of cases within each testing fold
    scaler = scaling_method()                                                   # defining the scaling method
    clf = SVC(kernel='linear', C=C)                                             # defining the classifier
    for t in range(X.shape[2]):
        x = X[:, :, t]                                                          # reducing data to a single time point
        case_inds = [i for i in range(m)]                                       # number of cases within each condition (assumed to be equal)
        for fold in range(k):
            test_inds = np.random.RandomState().choice(case_inds, fold_size,    # randomly sampling testing cases for this fold
                                                       False)
            case_inds = list(set(case_inds) - set(test_inds))                   # removing cases randomly sampled during line above from future pool of cases
            train_inds = np.array(list(set(np.arange(m)) - set(test_inds)))     # training cases are the complement of the testing cases
            assert len(test_inds)+len(train_inds) == m                          # making sure above logic is working correctly
            test_inds = np.hstack((test_inds, test_inds+m))                     # adding matched indicies for second cases (assumes that data is sorted so that i == i+m)
            train_inds = np.hstack((train_inds, train_inds+m))
            train_k = scaler.fit_transform(x[train_inds, :])                    # extracting fold's training data and fitting scaler
            train_y = y[train_inds]                                             # training labels
            test_k = scaler.transform(x[test_inds, :])                          # extraining fold's testing data and scaling
            test_y = y[test_inds]                                               # testing labels
            clf.fit(train_k, train_y)                                           # training the classifier
            predictions = clf.predict(test_k)                                   # testing the classifier
            accuracy = sum(predictions == test_y)/len(test_y)                   # calculating model accuracy
            output[t, fold] = accuracy                                          # storing classifier performance
            weights[t, :, fold] = project_weights(train_k, clf)

    return output.mean(axis=1), weights.mean(axis=2)

def linear_svm_TGM(X_1, X_2, y, m, C=1, k=10):
    """
    Function written for linear temporal generalisation decoding of evoked EEG 
    data using a support-vector classifier.
    
    Parameters
    ----------
    X_1 : TYPE
        DESCRIPTION.
    X_2 : TYPE
        DESCRIPTION.
    y : TYPE
        DESCRIPTION.
    m : TYPE
        DESCRIPTION.
    win : TYPE
        DESCRIPTION.
    step : TYPE
        DESCRIPTION.
    C : TYPE, optional
        DESCRIPTION. The default is 1.
    k : TYPE, optional
        DESCRIPTION. The default is 2.

    Returns
    -------
    TGM : array
        Temporal generalisation matrix data averaged over cross-validation 
        folds.
    """
    assert int(len(y)/len(list(set(y)))) == int(m), ('ERROR')
    assert X_1.shape == X_2.shape, ('ERROR')
    output = np.zeros((X_1.shape[2], X_1.shape[2], k))
    fold_size = floor(m/k)
    scaler = scaling_method()
    clf = SVC(kernel='linear', C=C)
    for t_1 in range(X_1.shape[2]):
        x_1 = X_1[:, :, t_1]
        for t_2 in range(X_2.shape[2]):
            x_2 = X_2[:, :, t_2]
            case_inds = [i for i in range(m)]
            assert x_1.shape == x_2.shape
            for fold in range(k):
                test_inds = np.random.RandomState().choice(case_inds, fold_size, False) # randomly determing fold's test inds
                case_inds = list(set(case_inds) - set(test_inds))               # removing above determined test inds from the pool of cases so that they can't be selected as testing indicies again during the next loop of k
                train_inds = np.array(list(set(np.arange(m)) - set(test_inds))) # determing training indicies as all of those that aren't testing indicies
                assert len(test_inds)+len(train_inds) == m, ('ERROR')           # asserting that the training and testing indicies add to m (the number of cases within each condition)
                test_inds = np.hstack((test_inds, test_inds+m))                 # determining testing indicies for matched pairs within second condition
                train_inds = np.hstack((train_inds, train_inds+m))              # determining training indicies for the matched pairs within the second condition                               # scaling data
                train_k = scaler.fit_transform(x_1[train_inds, :])              # scaling training data
                train_y = y[train_inds]                                         # training labels
                test_k = scaler.transform(x_2[test_inds, :])                    # scaling testing data by parameter fit two lines above
                test_y = y[test_inds]                                           # testomg labels
                clf.fit(train_k, train_y)                                       # fitting classifier
                predictions = clf.predict(test_k)                               # prediciting testing data
                accuracy = sum(predictions == test_y)/len(test_y)               # calculating accuracy
                output[t_1, t_2, fold] = accuracy                               # storing output

    return output.mean(axis=2)

#%% Decoding pipeline
def process_participant(fif_file, sub, conditions, baseline, ch_picks,
                        min_trials, behav_data, C=1, k=10, rep_avg=False):
    """
    Function written for running the analysis pipeline. This is written in a
    way that is somewhat akin to a main script (and probably should have gone
    in there) but within a function so that it could be called multiple times
    in parallel across CPU cores using ProcessPoolExecutor().

    See: https://docs.python.org/3/library/concurrent.futures.html
    """
    
    ## reading in data
    epochs = mne.read_epochs(fif_file, preload='True')
    epochs.apply_baseline(baseline)

    ## reducing data dimensionality to picked electrodes
    epochs = epochs.pick_channels(ch_picks, ordered=True)

    ## exactracting main effect data
    all_vio = epochs['col_vio', 'siz_vio', 'ori_vio'].get_data()
    control = epochs['control'].get_data()

    ## Removing unmatched epochs
    epochs.events = find_new_codes(epochs.events, behav_data)
    main_effect_data = get_main_effect_data(epochs, all_vio, control)
    del all_vio, control
    bad_epochs = find_unmatched_trials(epochs.events)
    epochs.drop(bad_epochs)

    ## Averaging trial repetitions
    if rep_avg:
        sorted_data = average_reps(epochs)
        m = sorted_data[conditions[0]].shape[0]
        assert m >= min_trials, "{} did not contain enough trials".format(sub)

    ## Reordering epoch data and sorting them into their correct condition (only implemented if repetitions were not averaged)
    else:
        m = int(len(epochs)/4)
        assert m >= min_trials, "{} did not contain enough trials".format(sub)
        new_epoch_order = find_reorder_indicies(epochs)
        data = epochs.get_data()
        data = data[new_epoch_order, :, :]
        sorted_data = {}
        for i, j in enumerate(range(0, data.shape[0], m)):                      # sorting ordered data into the correct condition
            inds = list(range(j, j+m))
            sorted_data[conditions[i]] = data[inds]

    ## Decoding violations from controls (TGM)
    vio_vs_cont = np.zeros((len(conditions)-1, len(epochs.times),
                            len(epochs.times)))
    for v, violation in enumerate(conditions[:3]):
        print(f'Classifying {violation} against control.')
        X = np.vstack((sorted_data[violation], sorted_data['control']))
        y = np.hstack((np.zeros(m), np.ones(m)))
        vio_vs_cont[v, :, :] = linear_svm_TGM(X, X, y, m, C, k)

    ## Weights analysis (violations vs. controls)
    vio_vs_cont_SL = np.zeros((2, len(epochs.times), len(epochs.ch_names)))
    for v, violation in enumerate(conditions[1:3]):
        print(f'Classifying {violation} against control (weights).')
        X = np.vstack((sorted_data[violation], sorted_data['control']))
        y = np.hstack((np.zeros(m), np.ones(m)))
        DELETE, vio_vs_cont_SL[v, :, :] = linear_svm_decode(X, y, m, C, k)

    ## Decoding size violations from orientation violations
    X = np.vstack((sorted_data['siz_vio'], sorted_data['ori_vio']))
    y = np.hstack((np.zeros(m), np.ones(m)))
    print('Classifying siz_vio against ori_vio.')
    siz_vs_ori = linear_svm_TGM(X, X, y, m, C, k)                               # TGM
    print('Classifying siz_vio against ori_vio (weights).')
    DELETE, siz_vs_ori_SL = linear_svm_decode(X, y, m, C, k)                    # diagonal weights

    ## Transfer learning
    TL = np.zeros((2, 2, len(epochs.times), len(epochs.times)))
    for v_1, violation_1 in enumerate(conditions[1:3]):
        X_train = np.vstack((sorted_data[violation_1], sorted_data['control']))
        y = np.hstack((np.zeros(m), np.ones(m)))
        for v_2, violation_2 in enumerate(conditions[1:3]):
            if v_1 == v_2:                                                      # ensuring generalised conditions are different from one another
                continue
            print(f'Training on {violation_1}, testing on {violation_2}.')
            X_test = np.vstack((sorted_data[violation_2], sorted_data['control']))
            TL[v_1, v_2, :, :] = linear_svm_TGM(X_train, X_test, y, m, C, k)

    ## Main effect of expectation (not yet properly matched across folds)
    main_effect = np.zeros((len(epochs.times), len(epochs.times)))
    main_effect_SL = np.zeros((len(epochs.times), len(epochs.ch_names)))
    X = np.vstack((main_effect_data['all_vio'], main_effect_data['control']))
    m = int(X.shape[0]/2)
    y = np.hstack((np.zeros(m), np.ones(m)))
    print('Classifying main effect')
    main_effect = linear_svm_TGM(X, X, y, m, C, k)
    DELETE, main_effect_SL = linear_svm_decode(X, y, m, C, k)

    return (vio_vs_cont, vio_vs_cont_SL, siz_vs_ori, siz_vs_ori_SL, TL,
            main_effect, main_effect_SL)
