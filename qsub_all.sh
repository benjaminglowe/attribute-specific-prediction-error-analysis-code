#!/bin/bash -I

## DESCRIPTION ##
# ============= #
# This was the only line of code run by the on QUT's HPC clusters following
# preprocessing.py. Effectively, it queues a decoding job (see qsub_SS.sh for
# job parameters) for each fif_file within the fif_files list.

# fif files scring
fif_files=('sub-01-epo.fif' 'sub-02-epo.fif' 'sub-03-epo.fif' 'sub-04-epo.fif'
           'sub-05-epo.fif' 'sub-06-epo.fif' 'sub-07-epo.fif' 'sub-08-epo.fif'
           'sub-10-epo.fif' 'sub-12-epo.fif' 'sub-13-epo.fif' 'sub-14-epo.fif'
           'sub-15-epo.fif' 'sub-17-epo.fif' 'sub-20-epo.fif' 'sub-21-epo.fif'
           'sub-22-epo.fif' 'sub-23-epo.fif' 'sub-24-epo.fif' 'sub-25-epo.fif'
           'sub-26-epo.fif' 'sub-28-epo.fif' 'sub-29-epo.fif' 'sub-30-epo.fif'
           'sub-32-epo.fif' 'sub-33-epo.fif' 'sub-34-epo.fif' 'sub-35-epo.fif'
           'sub-36-epo.fif' 'sub-37-epo.fif')

# qsub command
for fif_file in ${fif_files[*]}; do
    qsub -v input=$fif_file qsub_SS.sh
done
