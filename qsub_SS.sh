#!/bin/bash -l
#PBS -N EEG_decode_SS
#PBS -l ncpus=32
#PBS -l mem=64GB
#PBS -l walltime=168:00:00

## DESCRIPTION ##
# ============= #
# This is called by qsub_all.sh, and not directly by the user. Basically, it
# just queues a job on QUT's HPC clusters (collectively called 'Lyra') to run
# decoding_lyra.py.

## Determining fif file
fif_file=${input}
echo $fif_file

## Loading in python module
module load python/3.8.2-gcccore-9.3.0

## Resetting working directory
cd ~
cd datasets/EEG_decode_v2/code

## Running code
python decoding_lyra.py $fif_file

echo "code ran"
