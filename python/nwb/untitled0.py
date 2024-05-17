# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 16:05:50 2024

@author: zhixi
"""
#%%
import os
import numpy as np

import shutil
import json
from pathlib import Path

from uuid import uuid4

import probeinterface as pi
import spikeinterface as si
# needed to lead extensions
import spikeinterface.postprocessing as spost
import spikeinterface.qualitymetrics as sqm
  
from pynwb import NWBHDF5IO
from pynwb.file import Device
from hdmf_zarr import NWBZarrIO

# from utils import get_devices_from_metadata, add_waveforms_with_uneven_channels

#%%
session = '689514_2024-02-01_18-06-43'

# parse sessoin info
tmpStr = session.split('_')
baseFolder = 'F:\acuteBehavior'
data_folder = Path(os.path.join('F:\\acuteBehavior', tmpStr[0],'689514_2024-02-01_18-06-43'))

io_class = NWBHDF5IO


# find sorted data
sorted_folders = [
    p for p in data_folder.iterdir() if p.is_dir() and "sorted" in p.name and "spikesorted" not in p.name and session in p.name
    ]

# find nwb
nwb_files = [p for p in data_folder.glob("**/*") if p.name.endswith(".nwb")]
#%%
io = NWBHDF5IO(nwb_files[0], mode='r')
nwb = io.read()
#%%