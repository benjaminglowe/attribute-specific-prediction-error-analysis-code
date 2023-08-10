# Analysis code for 'Same but Different: The Latency of a Shared Expectation Signal Interacts with Stimulus Attributes'
Hi there and thank you for visiting this directory, which contains all code necessary to replicate our analysis (DOI LINK HERE).

We have done our best to ensure that each file is sufficiently commented so that users easily navigate our scripts, however, a general overview of what each file does is as follows:
- `EEG_channels_64.csv` simply contains some meta information concerning the names and locations of the EEG channels used during our experiment. We only ever call this during line 32 of `decoding_lyra.py`.
- `preprocessing.py` contains our EEG preprocessing code, which depends on the [MNE-Python library](https://mne.tools/stable/index.html) and `project_funcs.py`.
- `decoding_lyra.py` is our main script for MVPA analysis. It calls functions defined within `project_funcs.py`.
- `get_erp_data.py` outputs the event-related potential data used to generate Fig. 2--Supplement 3.

You will notice three bash scripts. These are included for the sake of transparency and are simply the files we used to queue jobs on QUT's High Performance Computer clusters (called "Lyra"):
- `qsub_preprocess_all.sh`: queued `preprocessing.py`.
- `qsub_SS.sh`: queued `decoding_lyra.py` for a single subject (hence, "SS")
- `qsub_all.sh`: wrapped `qsub_SS.sh` to queue `decoding_lyra.sh` for all specified subjects.

All Figures within the manuscript and supplementary material were generated using the MATLAB scripts within the `statistical_inference` directory.
- `results_violations_vs_controls.m`: Fig. 2.
- `results_weights.m`: Fig. 2--Supplment 2.
- `results_erp.m`: Fig. 2--Supplement 3.
- `results_siz_vs_ori.m`: Fig. 3.
- `results_transfer_learning_TGM.m`: Fig. 4A and Supplementary Analysis.
- `results_temporal_shift.m`: Fig. 4B.
- The `.mat` files contain parameters that are loaded in when running these scripts.
- `functions` is a directory containing MATLAB functions used during this phase of analysis.

Our recommended run order is:
1. `preprocessing.py`
2. `decoding_lyra.py` (set MVPA = True)
3. `decoding_lyra.py` (set MVPA = False)
4. `get_erp_data.py`
5. `statistical_inference/results_violations_vs_controls.m`
6. `statistical_inference/results_siz_vs_ori.m`
7. `statistical_inference/results_transfer_learning_TGM.m'
8. `statistical_inference/results_temporal_shift.m`
9. Then the rest as the user sees fit.
