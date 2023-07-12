#!/bin/bash -l
#PBS -N EEG_preprocess
#PBS -l ncpus=16
#PBS -l mem=32GB
#PBS -l walltime=3:00:00

## DESCRIPTION ##
# ============= #
# This code was run on QUT's HPC clusters before any decoding analysis. I did
# this so that I did not have to vsync preprocessed data onto my HPC drive.
# Basically, it just queues a job to run the preprocessing.py script.

## Loading in python module
module load python/3.9.5-gcccore-10.3.0

## Resetting working directory
cd ~
cd datasets/EEG_decode_v2/code

## Running code
python preprocessing.py
