# Analysis code for 'Same but Different: The Latency of a Shared Expectation Signal Interacts with Stimulus Attributes'
Hi there and thank you for visiting this directory, which contains all code necessary to replicate our analysis (DOI LINK HERE).

We have done our best to ensure that each file is sufficiently commented so that users easily navigate our scripts, however, a general overview of what each file does is as follows:
- `EEG_channels_64.csv` simply contains some meta information concerning the names and locations of the EEG channels used during our experiment. We only ever call this during line 32 of `decoding_lyra.py`.
- `preprocessing.py` contains our EEG preprocessing code, which depends on the [MNE-Python library](https://mne.tools/stable/index.html).
- 

Our recommended run order is:
