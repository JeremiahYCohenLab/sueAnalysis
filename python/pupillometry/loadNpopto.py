# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 21:49:55 2023

@author: zhixi
"""
#%%
import spikeinterface as si
import numpy as np
from pathlib import Path
from scipy.io import savemat
import os
from spikeinterface.curation import apply_sortingview_curation

#%%
sessionID = '684930_2023-09-28_12-44-15'
savePath = fr"F:/npOptoRecordings/{sessionID}/sorted/processed/";

samplestampsPath = os.path.join(savePath, 'samplestamps')
waveformsPath = os.path.join(savePath, 'waveforms')
sorting_asset = Path(f"F:/npOptoRecordings\{sessionID}\sorted\sorting_precurated\experiment1_Record Node 104#Neuropix-PXI-100.ProbeA-AP_recording1");
waveform_folder = Path(f"F:/npOptoRecordings\{sessionID}\sorted\postprocessed\experiment1_Record Node 104#Neuropix-PXI-100.ProbeA-AP_recording1");
sorting_output = si.load_extractor(sorting_asset);
we = si.WaveformExtractor.load_from_folder(waveform_folder, with_recording=False, sorting=sorting_output);
print(np.shape(sorting_output.unit_ids))
#%%
#extremum_channels = si.get_template_extremum_channel(we)

units_to_keep = sorting_output.unit_ids.flatten()
sorting_good = sorting_output.select_units(units_to_keep)
units = sorting_good.unit_ids
unit_ids = we.sorting.unit_ids

allWF = we.get_all_templates();

if not os.path.isdir(savePath):
    os.mkdir(savePath)

if not os.path.isdir(samplestampsPath):
    os.mkdir(samplestampsPath)
    
if not os.path.isdir(waveformsPath):
    os.mkdir(waveformsPath)  

for unit_ind, unit in enumerate(units):
    print(f'current unit {unit}')
    unit_waveform = we.get_waveforms(unit)
    print(unit_waveform.shape)
    sample_numbers = sorting_output.get_unit_spike_train(unit)    
    np.save(f'{waveformsPath}\waveforms_{unit}.npy', sample_numbers)              
    np.savetxt(f'{samplestampsPath}\samplestamps_{unit}.txt', sample_numbers)

#%%
mdic = {'waveforms', allWF}

savemat("F:\npOptoRecordings\{sessionID}\sorted\waveforms\waveforms.mat", mdic)

sample_numbers = sorting_output.get_unit_spike_train(unit)
unit_spike_times = timestamps[sample_numbers]


extremum_channels = si.get_template_extremum_channel(we)

#quality_metrics = pd.read_csv(f"/data/ecephys_{session}_sorted-ks2_5/postprocessed/experiment1_Record Node 104#Neuropix-PXI-100.ProbeA-AP_recording1/quality_metrics/metrics.csv")

units_to_keep = sorting_output.unit_ids.flatten()
print('Unit number all')
print(units_to_keep)

#units_to_keep = sorting_output.unit_ids[sorting_output.get_property("pass_qc")].flatten()
#units_to_keep = quality_metrics.query('isi_violations_ratio < 0.5 and amplitude_cutoff < 0.1 and presence_ratio > 0.95')
sorting_good = sorting_output.select_units(units_to_keep)

# prepare dataframe to store laser response data
laser_response_metrics = pd.DataFrame({'unit_id':units_to_keep})
noise_labels = np.where(laser_response_metrics['unit_id'].isin(noise_units), 'noise', 'good')
laser_response_metrics['noise_label'] = noise_labels

param_group = 'train5Hz'
trial_types = np.unique(event_ids.type)

# pre-adding all the columns that will contain lists so they can store objects
for ind, trial_type in enumerate(trial_types):
    laser_response_metrics.insert(2+5*ind,f'{trial_type}_{param_group}_all_pVals_unpaired','')
    laser_response_metrics.insert(3+5*ind,f'{trial_type}_{param_group}_all_pVals_paired','')
    laser_response_metrics.insert(4+5*ind,f'{trial_type}_{param_group}_all_latencies','')
    laser_response_metrics.insert(5+5*ind,f'{trial_type}_{param_group}_all_jitters','')
    laser_response_metrics.insert(6+5*ind,f'{trial_type}_{param_group}_all_reliability','')

laser_response_metrics = laser_response_metrics.astype(object)

print('Unit number')
print(np.shape(sorting_good.unit_ids))
#%%

from neuroconv.tools.spikeinterface.spikeinterface import (
    add_electrodes_info,
    add_units_table,
    get_electrode_group_indices
)
add_units_table


np.load('C:/Users/zhixi/Downloads/spikes.npy')