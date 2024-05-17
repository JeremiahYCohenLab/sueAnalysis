# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 15:36:18 2024

@author: zhixi
"""
#%%
# import numpy as np
# import pandas as pd

# _SECONDS_PER_TICK = 32e-6
# _payloadtypes = {
#                 1 : np.dtype(np.uint8),
#                 2 : np.dtype(np.uint16),
#                 4 : np.dtype(np.uint32),
#                 8 : np.dtype(np.uint64),
#                 129 : np.dtype(np.int8),
#                 130 : np.dtype(np.int16),
#                 132 : np.dtype(np.int32),
#                 136 : np.dtype(np.int64),
#                 68 : np.dtype(np.float32)
#                 }

# def read_harp_bin(file):

#     data = np.fromfile(file, dtype=np.uint8)

#     if len(data) == 0:
#         return None

#     stride = data[1] + 2
#     length = len(data) // stride
#     payloadsize = stride - 12
#     payloadtype = _payloadtypes[data[4] & ~0x10]
#     elementsize = payloadtype.itemsize
#     payloadshape = (length, payloadsize // elementsize)
#     seconds = np.ndarray(length, dtype=np.uint32, buffer=data, offset=5, strides=stride)
#     ticks = np.ndarray(length, dtype=np.uint16, buffer=data, offset=9, strides=stride)
#     seconds = ticks * _SECONDS_PER_TICK + seconds
#     payload = np.ndarray(
#         payloadshape,
#         dtype=payloadtype,
#         buffer=data, offset=11,
#         strides=(stride, elementsize))

#     if payload.shape[1] ==  1:
#         ret_pd = pd.DataFrame(payload, index=seconds, columns= ["Value"])
#         ret_pd.index.names = ['Seconds']

#     else:
#         ret_pd =  pd.DataFrame(payload, index=seconds)
#         ret_pd.index.names = ['Seconds']

#     return ret_pd
# %%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import seaborn as sns
import harp
import os
import csv
 
# %% parameters
root = r'Z:\ephys\temp\Behavior'
session = '0_2024-03-22_09-36-54';
videoTime = '09-37-39.130.ts.csv';
fileToRead = 'Event_94.bin';
# parse session name
cellTemp = session.split('_')
aniID = cellTemp[0]
date = cellTemp[1]
time = cellTemp[2]
file = os.path.join(root, aniID, session, 'HarpFolder', 'BehaviorEvents', fileToRead)
# metadata = pd.read_csv(root / r"VideoFolder\side_camera.csv", header=None)
triggers = harp.read(file)
# %% load frame time
csvFile = os.path.join(root, aniID, session, 'VideoFolder', videoTime)
a = pd.read_csv(csvFile)
#%%
import datetime
import time

FMT = '%H:%M:%S.%f'
timeRaw = a['+0:00:00.000']
tsTime = np.ones_like(timeRaw)
# calculate difference
for i in range(len(timeRaw)):
    currString = timeRaw[i][1:]
    sepTime = currString.split(':')
    tsTime[i] = 1000*(float(sepTime[0]) * 3600 + float(sepTime[1]) * 60 + float(sepTime[2]))
#%%
    
    
# output: 143.143
#%%