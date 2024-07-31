# -*- coding: utf-8 -*-

"""
Created on Tue Oct 10 21:49:55 2023

@author: zhixi
"""

import spikeinterface as si
import numpy as np
from pathlib import Path
from scipy.io import savemat

sessionID = '684930_2023-09-28_12-44-15'
savepath = f"F:\npOptoRecordings\{sessionID}\sorted\waveforms";
sorting_asset = Path(f"F:/npOptoRecordings\{sessionID}\sorted\spikesorted\experiment1_Record Node 104#Neuropix-PXI-100.ProbeA-AP_recording1");
waveform_folder = Path(f"F:/npOptoRecordings\{sessionID}\sorted\postprocessed\experiment1_Record Node 104#Neuropix-PXI-100.ProbeA-AP_recording1");
sorting_output = si.load_extractor(sorting_asset);
we = si.WaveformExtractor.load_from_folder(waveform_folder, with_recording=False, sorting=sorting_output);
print(np.shape(sorting_output.unit_ids))
print(len(sorting_output.unit_ids))
#%%
