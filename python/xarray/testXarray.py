# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 16:36:48 2024

@author: zhixi
"""
import xarray as xr
import numpy as np
import pandas as pd

import os
import matplotlib.pyplot as plt
from pynwb import NWBFile, TimeSeries, NWBHDF5IO

from pathlib import Path

#%% load allen data: 
# from allensdk.brain_observatory.ecephys.ecephys_project_cache import EcephysProjectCache
# from pathlib import Path
# manifest_path = r'F:\allenData'
# cache = EcephysProjectCache.from_warehouse(manifest=str(Path(manifest_path) / 'brain_observatory_manifest.json'))
# session = cache.get_session_data(819701982)  # access data by session ID
# #%% load stimulus
# stimulus_presentations = session.stimulus_presentations
# flash_presentations = stimulus_presentations[
#     stimulus_presentations.stimulus_name == "flashes"
# ]

# responses = session.presentationwise_spike_counts(
#     np.arange(0, 0.5, 0.001),
#     flash_presentations.index.values,
#     session.units.index.values
# )
# responses.coords

#%%

# nullArray = responses.values;
# time = np.arange(0, 0.5, 0.001)
# time = time[:-1]
# allResponses = xr.DataArray(nullArray, dims= ("trialInd", "time", "units"), 
#                             coords = {"trialInd" : np.array(range(150)), "time": time, "units": range(585)})
# desired_trial_inds = [1, 3, 5, 7]
# test = allResponses.sel(trialInd = desired_trial_inds, units = [2, 3, 4]).mean(dim='trialInd')
# test2 = allResponses.sel(trialInd = desired_trial_inds, units = [2, 3, 4]).mean(dim='trialInd').sortby("units").transpose()

#%% load nwb file
nwbFile = Path(r'F:\acuteBehavior\689514\689514_2024-02-01_18-06-43\sorted\689514_2024-02-01_18-06-43.nwb')
io = NWBHDF5IO(nwbFile, mode='r')
nwb = io.read()
tblTrials = nwb.trials.to_dataframe()
trialStarts = tblTrials['start_time'].values
responseTimes = tblTrials[tblTrials['animal_response']!=2]
responseTimes = responseTimes['reward_outcome_time'].values
#%% read example neuron timestamps
exampleFiles = [r"F:\acuteBehavior\689514\689514_2024-02-01_18-06-43\sorted\spiketimes\units_130.npy",
                r"F:\acuteBehavior\689514\689514_2024-02-01_18-06-43\sorted\spiketimes\units_132.npy",
                r"F:\acuteBehavior\689514\689514_2024-02-01_18-06-43\sorted\spiketimes\units_247.npy"]
timestampsRS = [np.load(a) for a in exampleFiles]
#%% all functions
def build_time_window_domain(bin_edges, offsets, callback=None):
    callback = (lambda x: x) if callback is None else callback
    domain = np.tile(bin_edges[None, :], (len(offsets), 1))
    domain += offsets[:, None]
    return callback(domain)

def build_spike_histogram(time_domain,
                          spike_times,
                          dtype=None,
                          binarize=False):

    time_domain = np.array(time_domain)

    tiled_data = np.zeros(
        (time_domain.shape[0], time_domain.shape[1] - 1, len(spike_times)),
        dtype=(np.uint8 if binarize else np.uint16) if dtype is None else dtype
    )

    starts = time_domain[:, :-1]
    ends = time_domain[:, 1:]
        
    for ii in range(len(spike_times)):
        data = np.array(spike_times[ii])

        start_positions = np.searchsorted(data, starts.flat)
        end_positions = np.searchsorted(data, ends.flat, side="right")
        counts = (end_positions - start_positions)

        tiled_data[:, :, ii].flat = counts > 0 if binarize else counts
    
    time = 0.5*(starts + ends)
    return time, tiled_data

def build_spike_histogram_overlap(time_domain,
                          binSize,
                          spike_times,
                          dtype=None,
                          binarize=False):

    time_domain = np.array(time_domain)

    tiled_data = np.zeros(
        (time_domain.shape[0], time_domain.shape[1] - 1, len(spike_times)),
        dtype=(np.uint8 if binarize else np.uint16) if dtype is None else dtype
    )

    starts = time_domain[:, :-1]
    ends = time_domain[:, :-1]  + binSize
        
    for ii in range(len(spike_times)):
        data = np.array(spike_times[ii])

        start_positions = np.searchsorted(data, starts.flat)
        end_positions = np.searchsorted(data, ends.flat, side="right")
        counts = (end_positions - start_positions)

        tiled_data[:, :, ii].flat = counts > 0 if binarize else counts

    time = starts + 0.5*binSize
    return time, tiled_data
#%% make ndarray for data array define functions first
bin_edges = np.arange(-2, 3, 0.1)
bin_edges = np.array(bin_edges)
time_domain_callback = None
alignTime = trialStarts
domain = build_time_window_domain(
    bin_edges,
    alignTime,
    callback=time_domain_callback)
dtype = None

ends = domain[:, -1]
starts = domain[:, 0]
time_diffs = starts[1:] - ends[:-1]
overlapping = np.where(time_diffs < 0)[0]

tiled_data = build_spike_histogram(
    domain,
    timestampsRS,
)

trialInds = np.array(range(len(alignTime)))

spike_toGoCue = xr.DataArray(
    name='spike_counts',
    data=tiled_data,
    dims=['trialInd',
          'trialTime',
          'unit_id'],
    coords={
        'trialInd': trialInds,
        'trialTime': (bin_edges[:-1] + np.diff(bin_edges) / 2),
        'unit_id': np.array(range(len(timestampsRS)))
    }
)
#%%









#%%

data_array = xr.DataArray(np.random.randn(4, 3), dims=('x', 'y'), coords={'x': [1, 2, 3, 4], 'y': ['a', 'b', 'c']})
print(data_array)

# Select specific rows and columns using isel
selected_data = data_array.isel(x=[0, 2], y=[1, 2])
print(selected_data)

#%%
